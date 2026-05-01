#!/usr/bin/env python3
from typing import Union

import h5py
import numpy as np


def combined_satellite_estimate(
    sub_pos,
    sub_dis,
    sub_index,
    sdss_direction,
    sdss_cos_surveyAngle,
    sdss_Rmax,
    des_direction,
    des_cos_surveyAngle,
    des_Rmax,
    MV,
    classical_sat_count,
    MV_bins,
):
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
            satellites could be there."""
            index_max = _index[select][0]
            N_tot[j] = np.random.randint(N_tot[j - 1], index_max + 1, 1)
            # This is the last iteration, so no need to get rid of any
            # subhaloes as the arrays won't be used again.

    # Add the classical satellites
    N_tot += classical_sat_count

    # Interpolate the satellite count on a regular grid in MV
    return 10 ** np.interp(MV_bins, MV, np.log10(N_tot))


def cartesian_to_spherical(coordinates: Union[list, np.ndarray]) -> np.ndarray:
    """Converts from Cartesian to Spherical basis.

    Args:
        coordinates (nd.array (N,3)): (x, y, z) values to convert.

    Returns:
        np.ndarray (N,3): coordinates in Spherical basis
            (r, theta, phi).
    """
    # Validation checks
    coordinates_permitted_types = (list, np.ndarray)
    if not isinstance(coordinates, coordinates_permitted_types):
        raise TypeError(
            "coordinates must be one of: {}".format(coordinates_permitted_types)
        )
    else:
        coordinates = np.asarray(coordinates)

    # Check if single vector is supplied or set of vectors
    array_of_arrays = isinstance(coordinates[0], np.ndarray)
    if not array_of_arrays:
        coordinates = coordinates.reshape((1, 3))

    # Function is designed to work in 3D
    if np.size(coordinates, 1) != 3:
        raise ValueError("coordinates should be (N,3) array")

    # Convert to spherical coordinates
    r = np.sqrt((coordinates * coordinates).sum(axis=1))
    theta = np.arccos(coordinates[:, 2] / r)
    phi = np.arctan2(coordinates[:, 1], coordinates[:, 0])

    # Output based on input array
    if array_of_arrays:
        return_array = np.column_stack((r, theta, phi))
    else:
        return_array = np.concatenate((r, theta, phi))

    return return_array


def compute200RhoCritical(M200, R200):
    """Returns the 200 * rho_critical value given M200 and R200.
            Courtesy of Marius Cautun.

    Args:
        M200 (fl) : Mass of halo [Msun/h].
        R200 (fl) : R_200 of halo [Mpc/h].

    Returns:
        fl: 200 * rho_critical.
    """
    return 3.0 / (4.0 * np.pi) * M200 / R200**3


def indiv_satellite_estimate(
    sub_pos,
    sub_dis,
    sub_index,
    direction,
    cos_surveyAngle,
    Rmax,
    MV,
    classical_sat_count,
    MV_bins,
):
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
            satellites could be there."""
            index_max = _index[select][0]
            N_tot[j] = np.random.randint(N_tot[j - 1], index_max + 1, 1)
            # This is the last iteration, so no need to get rid of any
            # subhaloes as the arrays won't be used again.

    # Add the classical satellites
    N_tot += classical_sat_count

    # Interpolate the satellite count on a regular grid in MV
    return 10 ** np.interp(MV_bins, MV, np.log10(N_tot))


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

    if cosTheta_min is None or cosTheta_min < -1.0:
        cosTheta_min = -1.0
    if cosTheta_max is None or cosTheta_max > +1.0:
        cosTheta_max = +1.0
    if cosTheta_min > cosTheta_max:
        cosTheta_min, cosTheta_max = cosTheta_max, cosTheta_min

    # z-coordinates
    sph[:, 2] = np.random.uniform(cosTheta_min, cosTheta_max, N)
    z2 = np.sqrt(1.0 - sph[:, 2] ** 2)
    phi = (2.0 * np.pi) * np.random.random(N)
    sph[:, 0] = z2 * np.cos(phi)  # x
    sph[:, 1] = z2 * np.sin(phi)  # y

    return sph


def second_pointings(pointings, alpha_off=2.0 * np.pi / 3.0):
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
    temp_w = w[w[:, 2] != 0.0]
    temp_u = u[w[:, 2] != 0.0]

    temp_u[:, 0] = 1.0
    temp_u[:, 1] = 1.0
    temp_u[:, 2] = -(temp_w[:, 0] + temp_w[:, 1]) / temp_w[:, 2]

    u[w[:, 2] != 0.0] = temp_u

    # wz == 0, wy != 0
    temp_w = w[(w[:, 2] == 0.0) & (w[:, 1] != 0.0)]
    temp_u = u[(w[:, 2] == 0.0) & (w[:, 1] != 0.0)]

    temp_u[:, 0] = 1.0
    temp_u[:, 1] = -(temp_w[:, 0] + temp_w[:, 2]) / temp_w[:, 1]
    temp_u[:, 2] = 1.0

    u[(w[:, 2] == 0.0) & (w[:, 1] != 0.0)] = temp_u

    # wz, wy ==0; wx !=0
    temp_w = w[(w[:, 1] == 0.0) & (w[:, 2] == 0.0)]
    temp_u = u[(w[:, 1] == 0.0) & (w[:, 2] == 0.0)]

    temp_u[:, 0] = -(temp_w[:, 1] + temp_w[:, 2]) / temp_w[:, 0]
    temp_u[:, 1] = 1.0
    temp_u[:, 2] = 1.0

    u[(w[:, 1] == 0.0) & (w[:, 2] == 0)] = temp_u

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
    new_pointings[:, 0] = (
        new_vecs[:, 0] * u[:, 0] + new_vecs[:, 1] * v[:, 0] + new_vecs[:, 2] * w[:, 0]
    )
    new_pointings[:, 1] = (
        new_vecs[:, 0] * u[:, 1] + new_vecs[:, 1] * v[:, 1] + new_vecs[:, 2] * w[:, 1]
    )
    new_pointings[:, 2] = (
        new_vecs[:, 0] * u[:, 2] + new_vecs[:, 1] * v[:, 2] + new_vecs[:, 2] * w[:, 2]
    )

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
    return 10 ** ((-a * Mv - b)) * 1.0e3


def survey_cone(survey_area):
    """Computes the fractional sky area coverage and cosine of a mock
        conical survey region.

    Args:
        survey_area (fl). Survey sky area coverage [deg^2].

    Returns:
        tuple:
            [0]: Fractional survey area
            [1]: Cosine(survey angle).
    """
    fullSkyArea = 4 * np.pi * (180.0 / np.pi) ** 2  # in deg^2
    fracArea = survey_area / fullSkyArea
    cos_surveyAngle = 1.0 - 2.0 * fracArea

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
    dz = 2.0 / N
    dPhi = np.pi * (3.0 - 5.0**0.5)
    points = np.empty((N, 3), np.float64)

    # Calculate z coordinate
    step = np.arange(N)
    points[:, 2] = 1.0 - dz * step

    # Calculate x and y coordinates
    r = (1.0 - points[:, 2] * points[:, 2]) ** 0.5
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
        out[i] = (MV_keep[i] >= MV_LMC).sum() + (MV_keep[i] >= classical_sats).sum()
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
    return (M200 / rhoCrit200 / (4.0 * np.pi / 3.0)) ** (1.0 / 3.0)


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
    with h5py.File(outputFile, "w") as hf:
        if groups:
            for g_i, g_item in enumerate(groups):
                grp = hf.create_group(g_item)
                grp.create_dataset("Lum_Func", data=data_out[g_i], dtype="f4")
        else:
            hf.create_dataset("Lum_Func", data=data_out[0], dtype="f4")
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
