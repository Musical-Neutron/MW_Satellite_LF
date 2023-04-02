#!/usr/bin/env python3
import sys
from typing import Union

import h5py
import numpy as np

survey_data = np.genfromtxt('Input_Data/Surveys.csv', names=True)
sat_data = np.genfromtxt('Input_Data/All_satellites.csv', names=True)
dic_surveyArea = {
    'sdss': survey_data['SDSSIXW09'][0],
    'des': survey_data['DES'][0]
}  # Survey sky area
# Column label in satellite file that gives SDSS and DES satellites
dic_columnSats = {'sdss': 'SDSSIX', 'des': 'DES'}

# a* and b* parameters to compute R_eff for a given magnitude.
dic_Reff_params = {
    'sdss': {
        'a': survey_data['SDSSIXW09'][1],
        'b': survey_data['SDSSIXW09'][2]
    },
    'des': {
        'a': survey_data['DES'][1],
        'b': survey_data['DES'][2]
    }
}
noSightings = 1000  # Default value for the number of sightings
M_max = -8.8  # Brightest magnitude bin used for interpolation grid
noBins_M = 45
M_min = M_max * (noBins_M / (0.5 - noBins_M) + 1.)
dM = M_max / (0.5 - noBins_M)
Sun_dis = 8.  # Distance of the sun from the Galactic Centre
disFactor = 1.e3  # Scale factor to obtain positions in 'kpc'
minVpeak = 10.  # Default value for the minimum V_peak used
maxRadius = 400.  # Default value for the maximum radius used
group_list = ['Original', 'Original+Orphan', 'Original+Orphan+Baryons']


def main():
    programName = sys.argv[0].rsplit('/')[-1]
    programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])
    help = '''
    ************************************************************************
    Computes the expected number of faint satellites given an observed
    number of satellites from a partial survey. This function requires 5
    command-line arguments in the below order. An additional 3 arguments
    may also be provided.

    1 : Input subhalo file  A post-processed subhalo file of the format
                            given in the documentation.
    2 : Host M200           The virial mass of the input host halo in
                            comoving coordinates (Msun h^-1).
    3 : Observation file    Observed satellites in the format given in the
                            documenation.
    4 : Target M200         The target (rescaled) virial mass of the host
                            halo in comoving coordinates (Msun h^-1).
    5 : Output file         The name of the output file to write luminosity
                            functions to.
    *6: N_sightings         The number of mock surveys per observer
                            position. DEFAULT: 1000.
    *7: Max radius          The fiducial radius inside which estimates are
                            made (in kpc). DEFAULT: 400.
    *8: Min vpeak           Imposes a cut greater than this value in subhalo
                            peak maximum circular velocity (in km/s).
                            DEFAULT: 10.

    where the * arguments are optional.
    NOTE: If optional arguments are given, they must be specified in the
    above order.
    ************************************************************************
    '''
    if len(sys.argv) not in [6, 7, 8, 9]:
        print(help)
        sys.exit(1)

    inputSubhaloFile = sys.argv[1]
    hostMass = float(sys.argv[2])
    inputObservFile = sys.argv[3]
    outputHostMass = float(sys.argv[4])
    outputFile = sys.argv[5]
    if len(sys.argv) >= 7:
        noSightings = int(sys.argv[6])
    if len(sys.argv) >= 8:
        maxRadius = float(sys.argv[7])
    if len(sys.argv) >= 9:
        minVpeak = float(sys.argv[8])
    """Computes the value of '200 rho_critical' needed to compute R200
    given a halo mass M200. The numerical values correspond to the M200
    and R200 values of Aq.A1."""
    delta200 = compute200RhoCritical(1.3432e12, 0.2458)

    # Obtain the fractional sky area/opening angle of each survey
    (sdss_surveyArea,
     sdss_cos_surveyAngle) = survey_cone(dic_surveyArea['sdss'])
    (des_surveyArea, des_cos_surveyAngle) = survey_cone(dic_surveyArea['des'])

    SDSS_data_out = np.empty(len(group_list)).tolist()
    DES_data_out = np.empty(len(group_list)).tolist()
    SDSSDES_data_out = np.empty(len(group_list)).tolist()

    for grp_i, grp_item in enumerate(group_list):
        # Read in the subhaloes
        print(
            "Reading subhalo data from file '{}' ...".format(inputSubhaloFile))
        with h5py.File(inputSubhaloFile, 'r') as data:
            pos = np.asarray(data[grp_item + '/cop'])
            vPeak = np.asarray(data[grp_item + '/vpeak'])

        # Select only objects above the Vpeak cut
        s = vPeak >= minVpeak
        s[0] = False  # Remove the host
        pos = pos[s] - pos[0]
        vPeak = vPeak[s]
        """ Compute a radial rescale factor. This is equivalent to changing
        the host mass to a new value. """
        # Input R200
        R200_original = virialRadius(hostMass, rhoCrit200=delta200)
        # R200 for desired output halo mass
        R200_output = virialRadius(outputHostMass, rhoCrit200=delta200)
        radialFactor = R200_output / R200_original
        print("Rescaling subhalo positions by a factor of {:0.2f} "
              "corresponding to changing from the input mass, M200={:0.2e}, "
              "to the desired output mass, M200 = {:0.2e}.".format(
                  radialFactor, hostMass, outputHostMass))

        # Select only subhaloes within the desired distance from the host
        dis = np.sqrt((pos * pos).sum(axis=1)) * disFactor * radialFactor
        # Select everything with distances <= maximum radius + Solar distance
        s = dis <= (maxRadius + Sun_dis)
        pos = pos[s, :] * disFactor * radialFactor
        vPeak = vPeak[s]

        # order the subhaloes according to their vMax values
        order = vPeak.argsort()[::-1]
        pos = pos[order, :]
        noSubs = pos.shape[0]
        print("There are {} subhaloes with Vmax >= {:0.1f} km/s and within a "
              "distance of {:0.0f} kpc from the central.".format(
                  noSubs, minVpeak, maxRadius + Sun_dis))
        """ Read the satellite data
        Reads magnitudes, distances, all provided satellite/survey data,
        and the probability of association with the LMC. """
        data = np.genfromtxt(inputObservFile, names=True)

        # SDSS satellites
        s = (data[dic_columnSats['sdss']] == 2) * (data['Dkpc'] <= maxRadius)
        MV_sdss = data['MV'][s]
        rad_sdss = data['Dkpc'][s]
        prob_sdss = data['LMC'][s]
        # Obtain classical satellites for SDSS
        s = (data[dic_columnSats['sdss']] == 1) * (data['Dkpc'] <= maxRadius)
        sdss_class = data['MV'][s]

        # DES satellites
        s = (data[dic_columnSats['des']] == 2) * (data['Dkpc'] <= maxRadius)
        MV_des = data['MV'][s]
        rad_des = data['Dkpc'][s]
        prob_des = data['LMC'][s]
        # Obtain classical satellites for DES
        s = (data[dic_columnSats['des']] == 1) * (data['Dkpc'] <= maxRadius)
        des_class = data['MV'][s]

        # Combine the two for the joint extrapolation
        MV_tot = np.append(MV_sdss, MV_des)
        rad_tot = np.append(rad_sdss, rad_des)
        prob_tot = np.append(prob_sdss, prob_des)
        # Obtain classical satellites for the combined SDSS and DES surveys
        s = ((data[dic_columnSats['sdss']] == 1) *
             (data['Dkpc'] <= maxRadius) + (data[dic_columnSats['des']] == 1) *
             (data['Dkpc'] <= maxRadius))
        tot_class = data['MV'][s]

        # Obtain 11 brightest classical satellites that are complete
        s = (data['Cl'] == 2) * (data['Dkpc'] <= maxRadius)
        class_class = data['MV'][s]

        print("There are {} satellites with magnitudes from {:0.1f} to "
              "{:0.1f}".format(MV_tot.shape[0], MV_tot.min(), MV_tot.max()))
        order = MV_sdss.argsort()
        MV_sdss = MV_sdss[order]
        prob_sdss = prob_sdss[order]

        order = MV_des.argsort()
        MV_des = MV_des[order]
        prob_des = prob_des[order]

        order = MV_tot.argsort()
        MV_tot = MV_tot[order]
        prob_tot = prob_tot[order]
        """ Generate a set of lines-of-sight (LOS) for the central direction
        of the survey, uniformly distributed over the surface of a sphere."""
        sdss_surveyDirs = uniformPointsOnSphereSurface(noSightings)
        des_surveyDirs = second_pointings(sdss_surveyDirs)

        # Compute the number of satellites using the 1st approach
        # Obtain magnitude bins
        bins_MV = np.linspace(M_max, M_min, noBins_M + 1, True)
        val_MV = 0.5 * (bins_MV[1:] + bins_MV[:-1])
        sdss_bin = np.zeros((6 * noSightings, noBins_M), np.float32)
        des_bin = np.zeros((6 * noSightings, noBins_M), np.float32)
        tot_bin = np.zeros((6 * noSightings, noBins_M), np.float32)

        # Obtain bins for the individual satellites
        sdss_MV = np.append(bins_MV[0] - 0.1,
                            np.append(MV_sdss, bins_MV[-1] + 0.1))
        sdss_Rmax = survey_Reff(sdss_MV,
                                a=dic_Reff_params['sdss']['a'],
                                b=dic_Reff_params['sdss']['b'],
                                surveyArea=sdss_surveyArea)
        sdss_prob = np.append(0., np.append(prob_sdss, 0.))

        des_MV = np.append(bins_MV[0], np.append(MV_des, bins_MV[-1] + 0.1))
        des_Rmax = survey_Reff(des_MV,
                               a=dic_Reff_params['des']['a'],
                               b=dic_Reff_params['des']['b'],
                               surveyArea=des_surveyArea)
        des_prob = np.append(0., np.append(prob_des, 0.))

        tot_MV = np.append(bins_MV[0], np.append(MV_tot, bins_MV[-1] + 0.1))
        tot_prob = np.append(0., np.append(prob_tot, 0.))
        tot_Rmax_sdss = survey_Reff(tot_MV,
                                    a=dic_Reff_params['sdss']['a'],
                                    b=dic_Reff_params['sdss']['b'],
                                    surveyArea=sdss_surveyArea)
        tot_Rmax_des = survey_Reff(tot_MV,
                                   a=dic_Reff_params['des']['a'],
                                   b=dic_Reff_params['des']['b'],
                                   surveyArea=des_surveyArea)

        print("Looping over 6 observer positions with {} sightings for "
              "each observer ...".format(noSightings))
        for k in range(6):  # Loop over observer positions
            # Obtain observer position
            pos_Sun = np.zeros(3, np.float32)
            """ For k=0, sets observer to x=+Sun_dis, k=1 -> x=-Sun_dis,
            k->2 y=+Sun_dis, and so on. """
            pos_Sun[int(k / 2)] = Sun_dis * (1 - 2 * (k % 2))
            print("Observer loop k = {} with Sun position "
                  "at {} ".format(k + 1, pos_Sun))

            # Translate subhalo positions relative to the observer position
            pos2 = pos + pos_Sun
            # Calculate distance of the subhaloes from the observer
            dis = np.sqrt((pos2 * pos2).sum(axis=1))
            select = dis <= maxRadius

            # select the subhaloes within the required radius
            pos2 = pos2[select]
            dis = dis[select]
            index = np.arange(dis.shape[0])

            for i in range(noSightings):  # Loop over the viewing directions
                # Randomises the subhalo list for each viewing direction.
                new_order = np.arange(len(pos2))
                np.random.shuffle(new_order)
                pos2 = pos2[new_order]
                dis = dis[new_order]
                """ SDSS """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(sdss_MV.shape[0]) >= sdss_prob
                """ If any satellites are associated to LMC, update the
            classical satellite count """
                class_sats = update_classical_satellites(
                    sdss_MV[select], sdss_MV[~select], sdss_class)
                sdss_bin[k * noSightings + i] = indiv_satellite_estimate(
                    pos2, dis, index, sdss_surveyDirs[i, :],
                    sdss_cos_surveyAngle, sdss_Rmax[select], sdss_MV[select],
                    class_sats, val_MV)
                """ DES """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(des_MV.shape[0]) >= des_prob
                """ If any satellites are associated to LMC, update the
                classical satellite count """
                class_sats = update_classical_satellites(
                    des_MV[select], des_MV[~select], des_class)
                des_bin[k * noSightings + i] = indiv_satellite_estimate(
                    pos2, dis, index, des_surveyDirs[i, :],
                    des_cos_surveyAngle, des_Rmax[select], des_MV[select],
                    class_sats, val_MV)
                """ SDSS + DES """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(tot_MV.shape[0]) >= tot_prob
                """ If any satellites are associated to LMC, update the
                classical satellite count """
                class_sats = update_classical_satellites(
                    tot_MV[select], tot_MV[~select], tot_class)
                tot_bin[k * noSightings + i] = combined_satellite_estimate(
                    pos2, dis, index, sdss_surveyDirs[i, :],
                    sdss_cos_surveyAngle, tot_Rmax_sdss[select],
                    des_surveyDirs[i, :], des_cos_surveyAngle,
                    tot_Rmax_des[select], tot_MV[select], class_sats, val_MV)

        # Compile results
        SDSS_data_out[grp_i] = output_data_format(val_MV, sdss_bin,
                                                  class_class)
        DES_data_out[grp_i] = output_data_format(val_MV, des_bin, class_class)
        SDSSDES_data_out[grp_i] = output_data_format(val_MV, tot_bin,
                                                     class_class)

    # SDSS
    writeOutput(outputFile + 'SDSS.hdf5', programOptionsDesc, SDSS_data_out,
                group_list)

    # DES
    writeOutput(outputFile + 'DES.hdf5', programOptionsDesc, DES_data_out,
                group_list)

    # DES + SDSS
    writeOutput(outputFile + 'SDSS+DES.hdf5', programOptionsDesc,
                SDSSDES_data_out, group_list)

    return None


def combined_satellite_estimate(sub_pos, sub_dis, sub_index, sdss_direction,
                                sdss_cos_surveyAngle, sdss_Rmax, des_direction,
                                des_cos_surveyAngle, des_Rmax, MV,
                                classical_sat_count, MV_bins):
    """Creates a luminosity function estimate for a combined observation
        of two mock surveys.

    Args:
        sub_pos (Nx3 arr): Cartesian subhalo positions [kpc].
        sub_dis (Nx1 arr): Subhalo distances from halo centre [kpc].
        sub_index (Nx1 arr): Indexes that order the subhalo list.
        sdss_direction (1x3 arr): Vector specifying the direction of the
            SDSS mock survey cone.
        sdss_cos_surveyAngle (fl): Cosine of the SDSS survey opening
            angle.
        sdss_Rmax (1xM arr): The SDSS Rmax value calculated for the
            associated M_V.
        des_direction (1x3 arr): Vector specifying the direction of the
            DES mock survey cone.
        des_cos_surveyAngle (fl): Cosine of the DES survey opening
            angle.
        des_Rmax (1xM arr): The DES Rmax value calculated for the
            associated M_V.
        MV (1xM arr): Absolute V-band magnitudes to be used in the
            estimate.
        classical_sat_count: Output from 'update classical satellites'.
        MV_bins (arr): Magnitude bins for estimate.

    Returns:
        arr: Satellite galaxy luminosity function.
    """
    # Find all subhaloes inside the mock surveys
    sdss_cos_theta = (sub_pos * sdss_direction).sum(axis=1) / sub_dis
    des_cos_theta = (sub_pos * des_direction).sum(axis=1) / sub_dis
    # Subhaloes inside the SDSS mock survey volume
    s_sdss = sdss_cos_theta >= sdss_cos_surveyAngle
    # Subhaloes inside the DES mock survey volume
    s_des = des_cos_theta >= des_cos_surveyAngle
    # Select subhaloes that are in either the first or second mocks
    select = s_sdss + s_des

    # Construct arrays of distances and indices for all subhaloes inside
    # the survey volumes.
    _dis = sub_dis[select]
    _index = sub_index[select]
    _is_sdss = s_sdss[select]  # True if subhalo is inside the SDSS mock

    # Temporary array to store the total number of subhaloes at each
    # observed satellite magnitude.
    N_tot = np.zeros(MV.shape[0], np.float32)
    """ Loop over each observed satellite. The first and last entries of
    the array are not actual observations so skip them. """
    for j in range(1, MV.shape[0]):
        # Select subhaloes inside SDSS and DES, respectively.
        s_sdss = +_is_sdss * (_dis <= sdss_Rmax[j])
        s_des = ~_is_sdss * (_dis <= des_Rmax[j])
        select = s_sdss + s_des
        if select.sum() < 2:
            # Check that we have at least two subhaloes for given satellite
            # (except for the last value of j)
            N_tot[j] = np.nan
        elif j < MV.shape[0] - 1:
            # Typical case, where we have 1 observation
            index_min, index_max = _index[select][[0, 1]]
            N_tot[j] = np.random.randint(index_min + 1, index_max + 1, 1)

            # Only keep subhaloes NOT used for the current estimate
            s2 = _index >= N_tot[j]
            # Update the _dis and _index array to reflect that objects were
            # removed.
            _dis = _dis[s2]
            _index = _index[s2]
            _is_sdss = _is_sdss[s2]
        else:
            """The last value of j, corresponding to the left-most point on
            the graph. No observation so we make an estimate of how many
            satellites could be there. """
            index_max = _index[select][0]
            N_tot[j] = np.random.randint(N_tot[j - 1], index_max + 1, 1)
            # This is the last iteration, so no need to get rid of any
            # subhaloes as the arrays won't be used again.

    # Add the classical satellites
    N_tot += classical_sat_count

    # Interpolate the satellite count on a regular grid in MV
    return 10**np.interp(MV_bins, MV, np.log10(N_tot))


def spherical_to_cartesian(coordinates: Union[list, np.ndarray]) -> np.ndarray:
    """Converts from Spherical to Cartesian basis.

        Args:
            coordinates (nd.array (N,3)): (r, theta, phi) values to
                convert.

        Returns:
            np.ndarray (N,3): coordinates in Cartesian basis (x, y, z).
        """
    # Validation checks
    coordinates_permitted_types = (list, np.ndarray)
    if not isinstance(coordinates, coordinates_permitted_types):
        raise TypeError("coordinates must be one of: {}".format(
            coordinates_permitted_types))
    else:
        coordinates = np.asarray(coordinates)

    # Check if single vector is supplied or set of vectors
    array_of_arrays = isinstance(coordinates[0], np.ndarray)
    if not array_of_arrays:
        coordinates = coordinates.reshape((1, 3))

    # Function is designed to work in 3D
    if np.size(coordinates, 1) != 3:
        raise ValueError("coordinates should be (N,3) array")

    # Convert to Cartesian coordinates
    x = coordinates[:, 0] * np.sin(coordinates[:, 1]) * np.cos(coordinates[:,
                                                                           2])
    y = coordinates[:, 0] * np.sin(coordinates[:, 1]) * np.sin(coordinates[:,
                                                                           2])
    z = coordinates[:, 0] * np.cos(coordinates[:, 1])

    # Output based on input array
    if array_of_arrays:
        return_array = np.column_stack((x, y, z))
    else:
        return_array = np.concatenate((x, y, z))

    return return_array


def compute200RhoCritical(M200, R200):
    """ Returns the 200 * rho_critical value given M200 and R200.
            Courtesy of Marius Cautun.

    Args:
        M200 (fl) : Mass of halo [Msun/h].
        R200 (fl) : R_200 of halo [Mpc/h].

    Returns:
        fl: 200 * rho_critical.
    """
    return 3. / (4. * np.pi) * M200 / R200**3


def indiv_satellite_estimate(sub_pos, sub_dis, sub_index, direction,
                             cos_surveyAngle, Rmax, MV, classical_sat_count,
                             MV_bins):
    """Creates a luminosity function estimate for a single mock survey.

    Args:
        sub_pos (Nx3 arr): Cartesian subhalo positions [kpc].
        sub_dis (Nx1 arr): Subhalo distances from halo centre [kpc].
        sub_index (Nx1 arr): Indexes that order the subhalo list.
        direction (1x3 arr): Vector specifying the direction of the mock
            survey cone.
        cos_surveyAngle (fl): Cosine of the survey opening angle.
        Rmax (1xM arr): The Rmax value calculated for the associated
            M_V.
        MV (1xM arr): Absolute V-band magnitudes to be used in the
            estimate.
        classical_sat_count: Output from 'update classical satellites'.
        MV_bins (arr): Magnitude bins for estimate.

    Returns:
        arr: Satellite galaxy luminosity function.
    """
    # Find all subhaloes inside the mock survey
    cos_theta = (sub_pos * direction).sum(axis=1) / sub_dis
    select = cos_theta >= cos_surveyAngle

    # Construct arrays of distances and indices for all subhaloes inside
    # the survey.
    _dis = sub_dis[select]
    _index = sub_index[select]

    # Temporary array to store the total number of subhaloes at each
    # observed satellite magnitude.
    N_tot = np.zeros(MV.shape[0], np.float32)
    """ Loop over each observed satellite. The first and last entries of
    the array are not actual observations so skip them. """
    for j in range(1, MV.shape[0]):
        # Select subhaloes inside survey volume.
        select = _dis <= Rmax[j]
        if select.sum() < 2:
            # Check that we have at least two subhaloes for given satellite
            # (except for the last value of j)
            N_tot[j] = np.nan
        elif j < MV.shape[0] - 1:
            # Typical case, where we have 1 observation
            index_min, index_max = _index[select][[0, 1]]
            N_tot[j] = np.random.randint(index_min + 1, index_max + 1, 1)

            # Only keep subhaloes NOT used for the current estimate
            s2 = _index >= N_tot[j]
            # Update the _dis and _index array to reflect that objects were
            # removed.
            _dis = _dis[s2]
            _index = _index[s2]
        else:
            """The last value of j, corresponding to the left-most point on
            the graph. No observation so we make an estimate of how many
            satellites could be there. """
            index_max = _index[select][0]
            N_tot[j] = np.random.randint(N_tot[j - 1], index_max + 1, 1)
            # This is the last iteration, so no need to get rid of any
            # subhaloes as the arrays won't be used again.

    # Add the classical satellites
    N_tot += classical_sat_count

    # Interpolate the satellite count on a regular grid in MV
    return 10**np.interp(MV_bins, MV, np.log10(N_tot))


def output_data_format(MV_bins, N_bins, class_class):
    """Prepares the data for writing.

    Args:
        MV_bins (arr): Magnitudes bins used in the estimation.
        N_bins (arr): Partially-complete array of estimated+classical
            luminosity functions.
        class_class (arr): Magnitudes of the classical satellites to
            prepend.

    Returns:
        arr (N_bins + 1): Luminosity functions. The first column
            corresponds to the absolute V-band magnitude values.
    """
    noRows = MV_bins.shape[0] + class_class.shape[0]
    noColumns = 1 + N_bins.shape[0]
    out = np.empty((noRows, noColumns), np.float32)

    # Insert the classical satellites
    class_class.sort()
    noClass = class_class.shape[0]
    # Magnitude (Mv) values corresponding to the classical satellites
    out[:noClass, 0] = class_class
    # Prepend cumulative classical satellite count
    out[:noClass, 1:] = (np.arange(noClass) + 1)[:, np.newaxis]

    # Append the estimated satellite count
    out[noClass:, 0] = MV_bins  # Magnitude bins
    out[noClass:, 1:] = N_bins.transpose()  # Estimated satellite count

    return out


def randomPointsOnSphereSurface(N, cosTheta_min=None, cosTheta_max=None):
    """Generates an array of random points on the surface of a unit
    sphere.

    Args:
        N (int): Number of points to generate.
        cosTheta_min (fl, optional): Lower bound on cosTheta values.
            Defaults to None.
        cosTheta_max (fl, optional): Upper bound on cosTheta values.
            Defaults to None.

    Returns:
        Nx3 arr: Cartesian vector for each point.
    """
    sph = np.empty((N, 3), np.float)

    if cosTheta_min is None or cosTheta_min < -1.:
        cosTheta_min = -1.
    if cosTheta_max is None or cosTheta_max > +1.:
        cosTheta_max = +1.
    if cosTheta_min > cosTheta_max:
        cosTheta_min, cosTheta_max = cosTheta_max, cosTheta_min

    # z-coordinates
    sph[:, 2] = np.random.uniform(cosTheta_min, cosTheta_max, N)
    z2 = np.sqrt(1.0 - sph[:, 2]**2)
    phi = (2.0 * np.pi) * np.random.random(N)
    sph[:, 0] = z2 * np.cos(phi)  # x
    sph[:, 1] = z2 * np.sin(phi)  # y

    return sph


def second_pointings(pointings, alpha_off=2. * np.pi / 3.):
    """Generates a set of pointing vectors with respect to an initial
        set.

    Args:
        pointings (Nx3 arr): Cartesian vectors with respect to which to
            generate the second set of vectors.
        alpha_off (fl, optional): Offset angle (in radians) with respect
            to the intial vectors provided to the function.
            Defaults to 2.*np.pi/3 (120 degrees).

    Returns:
        Nx3 arr: Cartesian vectors.
    """

    # Define new (u, v, w) orthogonal basis, with w aligned along the
    # initial pointing directions (pointings).
    # Define random rotations around w axis
    new_vec_orientation = np.random.rand(len(pointings)) * 2 * np.pi
    new_vecs = np.empty((len(pointings), 3))

    # Populate new_vecs in transformed (u, v, w) basis
    new_vecs[:, 0] = np.sin(new_vec_orientation) * np.sin(alpha_off)
    new_vecs[:, 1] = -np.cos(new_vec_orientation) * np.sin(alpha_off)
    new_vecs[:, 2] = np.cos(alpha_off)

    # Compute transformation matrix to cartesian basis
    mod_pointings = np.linalg.norm(pointings, axis=1)

    w = pointings / mod_pointings[:, np.newaxis]

    u = np.empty((len(pointings), 3), np.float64)
    v = np.empty((len(pointings), 3), np.float64)

    # wz != 0
    temp_w = w[w[:, 2] != 0.]
    temp_u = u[w[:, 2] != 0.]

    temp_u[:, 0] = 1.
    temp_u[:, 1] = 1.
    temp_u[:, 2] = -(temp_w[:, 0] + temp_w[:, 1]) / temp_w[:, 2]

    u[w[:, 2] != 0.] = temp_u

    # wz == 0, wy != 0
    temp_w = w[(w[:, 2] == 0.) & (w[:, 1] != 0.)]
    temp_u = u[(w[:, 2] == 0.) & (w[:, 1] != 0.)]

    temp_u[:, 0] = 1.
    temp_u[:, 1] = -(temp_w[:, 0] + temp_w[:, 2]) / temp_w[:, 1]
    temp_u[:, 2] = 1.

    u[(w[:, 2] == 0.) & (w[:, 1] != 0.)] = temp_u

    # wz, wy ==0; wx !=0
    temp_w = w[(w[:, 1] == 0.) & (w[:, 2] == 0.)]
    temp_u = u[(w[:, 1] == 0.) & (w[:, 2] == 0.)]

    temp_u[:, 0] = -(temp_w[:, 1] + temp_w[:, 2]) / temp_w[:, 0]
    temp_u[:, 1] = 1.
    temp_u[:, 2] = 1.

    u[(w[:, 1] == 0.) & (w[:, 2] == 0)] = temp_u

    v[:, 0] = u[:, 1] * w[:, 2] - w[:, 1] * u[:, 2]
    v[:, 1] = u[:, 2] * w[:, 0] - w[:, 2] * u[:, 0]
    v[:, 2] = u[:, 0] * w[:, 1] - w[:, 0] * u[:, 1]

    # Normalise vectors
    mod_u = np.linalg.norm(u, axis=1)
    u /= mod_u[:, np.newaxis]
    mod_v = np.linalg.norm(v, axis=1)
    v /= mod_v[:, np.newaxis]

    # Compute new pointing vectors in Cartesian basis
    new_pointings = np.empty((len(pointings), 3), np.float64)
    new_pointings[:, 0] = (new_vecs[:, 0] * u[:, 0] +
                           new_vecs[:, 1] * v[:, 0] + new_vecs[:, 2] * w[:, 0])
    new_pointings[:, 1] = (new_vecs[:, 0] * u[:, 1] +
                           new_vecs[:, 1] * v[:, 1] + new_vecs[:, 2] * w[:, 1])
    new_pointings[:, 2] = (new_vecs[:, 0] * u[:, 2] +
                           new_vecs[:, 1] * v[:, 2] + new_vecs[:, 2] * w[:, 2])

    # Compute unit vectors of new pointings
    mod_new_pointings = np.linalg.norm(new_pointings, axis=1)
    new_pointings /= mod_new_pointings[:, np.newaxis]

    # Compute original pointing unit vectors
    pointings /= mod_pointings[:, np.newaxis]

    return new_pointings


def survey_Reff(Mv, a, b, surveyArea):
    """Computes the effective radius of detection of satellites of a
        given magnitude in a given survey.

    Args:
        Mv (arr): Absolute V-band mangitude values.
        a (fl): a* parameter from eq. 3 in Newton et al. (2018)
        b (fl): b* parameter from eq. 3 in Newton et al. (2018)
        surveyArea (fl): Fractional sky area covered by given survey.

    Returns:
        arr: Reff [kpc].
    """
    return 10**((-a * Mv - b)) * 1.e3


def survey_cone(survey_area):
    """ Computes the fractional sky area coverage and cosine of a mock
        conical survey region.

    Args:
        survey_area (fl). Survey sky area coverage [deg^2].

    Returns:
        tuple:
            [0]: Fractional survey area
            [1]: Cosine(survey angle).
    """
    fullSkyArea = 4 * np.pi * (180. / np.pi)**2  # in deg^2
    fracArea = survey_area / fullSkyArea
    cos_surveyAngle = 1. - 2. * fracArea

    return (fracArea, cos_surveyAngle)


def uniformPointsOnSphereSurface(N):
    """Generates an array of points uniformly distributed on the surface
        of a unit sphere.
        http://bit.ly/2cQClOc
        Courtesy of Marius Cautun.

    Args:
        N (int): Number of points to generate.

    Returns:
        Nx3 arr: Cartesian vectors for the generated points.
    """
    dz = 2. / N
    dPhi = np.pi * (3. - 5.**0.5)
    points = np.empty((N, 3), np.float64)

    # Obtain z coordinate
    step = np.arange(N)
    points[:, 2] = 1. - dz * step
    # Obtain x and y coordinates
    r = (1. - points[:, 2] * points[:, 2])**0.5
    phi = dPhi * step
    points[:, 0] = r * np.cos(phi)
    points[:, 1] = r * np.sin(phi)

    return points


def update_classical_satellites(MV_keep, MV_LMC, classical_sats):
    """Updates the count of classical satellites given the selection of
        satellites that are associated to the LMC.

    Args:
        MV_keep (arr): Absolute V-band mangnitudes of the satellites to
            use in the procedure.
        MV_LMC (arr): Absolute V-band magnitudes of the satellites to be
            removed from the procedure.
        classical_sats (arr): Classical satellites to add in.

    Returns:
        1xN arr: Cumulative number of satellites that have been treated
        'classically' in this analysis.
    """
    out = np.zeros(MV_keep.shape[0], np.int32)
    for i in range(MV_keep.shape[0]):
        out[i] = ((MV_keep[i] >= MV_LMC).sum() +
                  (MV_keep[i] >= classical_sats).sum())
    return out


def virialRadius(M200, rhoCrit200):
    """Calculates R_200 of a halo given M_200 and rho_crit.

    Args:
        M200 (fl): Mass of halo [Msun / h].
        rhoCrit200 (fl): Critical density of the Universe
            [Msun/h (Mpc/h)^-3].

    Returns:
        fl: R_200 [Mpc / h].
    """
    return (M200 / rhoCrit200 / (4. * np.pi / 3.))**(1. / 3.)


def writeOutput(outputFile, programOptionsDesc, data_out, groups=False):
    """Writes 'Lum_Func' data set to hdf5 file. This has the format
        returned by 'output_data_format'.

    Args:
        outputFile (str): Full path to output file.
        programOptionsDesc (str): Command line command used to produce
            this set of luminosity functions.
        data_out (tuple): The data to write to the file.
        groups (bool, optional): If list, list of groups to write into
            output file. Defaults to False.

    Returns:
        None
    """
    print("Writing the output data file '{}' ...".format(outputFile))
    with h5py.File(outputFile, 'w') as hf:
        if groups:
            for g_i, g_item in enumerate(groups):
                grp = hf.create_group(g_item)
                grp.create_dataset('Lum_Func', data=data_out[g_i], dtype='f4')
        else:
            hf.create_dataset('Lum_Func', data=data_out[0], dtype='f4')
        grp = hf.create_group("Metadata")
        grp.attrs["Program options"] = np.string_(programOptionsDesc)
        nOptions = programOptionsDesc.split()
        if len(nOptions) >= 7:
            grp.attrs["Number of sightings"] = nOptions[6]
        else:
            grp.attrs["Number of sightings"] = noSightings
        if len(nOptions) >= 8:
            grp.attrs["Maximum radius"] = nOptions[7]
        else:
            grp.attrs["Maximum radius"] = maxRadius
        if len(nOptions) >= 9:
            grp.attrs["Minimum V_peak"] = nOptions[8]
        else:
            grp.attrs["Minimum V_peak"] = minVpeak
    return None


# Run script
if __name__ == "__main__":
    main()
