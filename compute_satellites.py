#!/usr/bin/env python3
import sys

import h5py
import numpy as np
from mwsatslf.tools import (
    combined_satellite_estimate,
    compute200RhoCritical,
    indiv_satellite_estimate,
    output_data_format,
    second_pointings,
    survey_Reff,
    survey_cone,
    uniformPointsOnSphereSurface,
    update_classical_satellites,
    virialRadius,
    writeOutput,
)

survey_data = np.genfromtxt("Input_Data/Surveys.csv", names=True)
sat_data = np.genfromtxt("Input_Data/All_satellites.csv", names=True)
dic_surveyArea = {
    "sdss": survey_data["SDSSIXW09"][0],
    "des": survey_data["DES"][0],
}  # Survey sky area
# Column label in satellite file that gives SDSS and DES satellites
dic_columnSats = {"sdss": "SDSSIX", "des": "DES"}

# a* and b* parameters to compute R_eff for a given magnitude.
dic_Reff_params = {
    "sdss": {"a": survey_data["SDSSIXW09"][1], "b": survey_data["SDSSIXW09"][2]},
    "des": {"a": survey_data["DES"][1], "b": survey_data["DES"][2]},
}
noSightings = 1000  # Default value for the number of sightings
M_max = -8.8  # Brightest magnitude bin used for interpolation grid
noBins_M = 45
M_min = M_max * (noBins_M / (0.5 - noBins_M) + 1.0)
dM = M_max / (0.5 - noBins_M)
Sun_dis = 8.0  # Distance of the sun from the Galactic Centre
disFactor = 1.0e3  # Scale factor to obtain positions in 'kpc'
minVpeak = 10.0  # Default value for the minimum V_peak used
maxRadius = 400.0  # Default value for the maximum radius used
group_list = ["Original", "Original+Orphan", "Original+Orphan+Baryons"]


def main():
    programName = sys.argv[0].rsplit("/")[-1]
    programOptionsDesc = programName + " " + " ".join(sys.argv[1:])
    help = """
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
    """
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
    sdss_surveyArea, sdss_cos_surveyAngle = survey_cone(dic_surveyArea["sdss"])
    des_surveyArea, des_cos_surveyAngle = survey_cone(dic_surveyArea["des"])

    SDSS_data_out = np.empty(len(group_list)).tolist()
    DES_data_out = np.empty(len(group_list)).tolist()
    SDSSDES_data_out = np.empty(len(group_list)).tolist()

    for grp_i, grp_item in enumerate(group_list):
        # Read in the subhaloes
        print("Reading subhalo data from file '{}' ...".format(inputSubhaloFile))
        with h5py.File(inputSubhaloFile, "r") as data:
            pos = np.asarray(data[grp_item + "/cop"])
            vPeak = np.asarray(data[grp_item + "/vpeak"])

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
        print(
            "Rescaling subhalo positions by a factor of {:0.2f} "
            "corresponding to changing from the input mass, M200={:0.2e}, "
            "to the desired output mass, M200 = {:0.2e}.".format(
                radialFactor, hostMass, outputHostMass
            )
        )

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
        print(
            "There are {} subhaloes with Vmax >= {:0.1f} km/s and within a "
            "distance of {:0.0f} kpc from the central.".format(
                noSubs, minVpeak, maxRadius + Sun_dis
            )
        )
        """ Read the satellite data
        Reads magnitudes, distances, all provided satellite/survey data,
        and the probability of association with the LMC. """
        data = np.genfromtxt(inputObservFile, names=True)

        # SDSS satellites
        s = (data[dic_columnSats["sdss"]] == 2) * (data["Dkpc"] <= maxRadius)
        MV_sdss = data["MV"][s]
        rad_sdss = data["Dkpc"][s]
        prob_sdss = data["LMC"][s]
        # Obtain classical satellites for SDSS
        s = (data[dic_columnSats["sdss"]] == 1) * (data["Dkpc"] <= maxRadius)
        sdss_class = data["MV"][s]

        # DES satellites
        s = (data[dic_columnSats["des"]] == 2) * (data["Dkpc"] <= maxRadius)
        MV_des = data["MV"][s]
        rad_des = data["Dkpc"][s]
        prob_des = data["LMC"][s]
        # Obtain classical satellites for DES
        s = (data[dic_columnSats["des"]] == 1) * (data["Dkpc"] <= maxRadius)
        des_class = data["MV"][s]

        # Combine the two for the joint extrapolation
        MV_tot = np.append(MV_sdss, MV_des)
        rad_tot = np.append(rad_sdss, rad_des)
        prob_tot = np.append(prob_sdss, prob_des)
        # Obtain classical satellites for the combined SDSS and DES surveys
        s = (data[dic_columnSats["sdss"]] == 1) * (data["Dkpc"] <= maxRadius) + (
            data[dic_columnSats["des"]] == 1
        ) * (data["Dkpc"] <= maxRadius)
        tot_class = data["MV"][s]

        # Obtain 11 brightest classical satellites that are complete
        s = (data["Cl"] == 2) * (data["Dkpc"] <= maxRadius)
        class_class = data["MV"][s]

        print(
            "There are {} satellites with magnitudes from {:0.1f} to "
            "{:0.1f}".format(MV_tot.shape[0], MV_tot.min(), MV_tot.max())
        )
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
        sdss_MV = np.append(bins_MV[0] - 0.1, np.append(MV_sdss, bins_MV[-1] + 0.1))
        sdss_Rmax = survey_Reff(
            sdss_MV,
            a=dic_Reff_params["sdss"]["a"],
            b=dic_Reff_params["sdss"]["b"],
            surveyArea=sdss_surveyArea,
        )
        sdss_prob = np.append(0.0, np.append(prob_sdss, 0.0))

        des_MV = np.append(bins_MV[0], np.append(MV_des, bins_MV[-1] + 0.1))
        des_Rmax = survey_Reff(
            des_MV,
            a=dic_Reff_params["des"]["a"],
            b=dic_Reff_params["des"]["b"],
            surveyArea=des_surveyArea,
        )
        des_prob = np.append(0.0, np.append(prob_des, 0.0))

        tot_MV = np.append(bins_MV[0], np.append(MV_tot, bins_MV[-1] + 0.1))
        tot_prob = np.append(0.0, np.append(prob_tot, 0.0))
        tot_Rmax_sdss = survey_Reff(
            tot_MV,
            a=dic_Reff_params["sdss"]["a"],
            b=dic_Reff_params["sdss"]["b"],
            surveyArea=sdss_surveyArea,
        )
        tot_Rmax_des = survey_Reff(
            tot_MV,
            a=dic_Reff_params["des"]["a"],
            b=dic_Reff_params["des"]["b"],
            surveyArea=des_surveyArea,
        )

        print(
            "Looping over 6 observer positions with {} sightings for "
            "each observer ...".format(noSightings)
        )
        for k in range(6):  # Loop over observer positions
            # Obtain observer position
            pos_Sun = np.zeros(3, np.float32)
            """ For k=0, sets observer to x=+Sun_dis, k=1 -> x=-Sun_dis,
            k->2 y=+Sun_dis, and so on. """
            pos_Sun[int(k / 2)] = Sun_dis * (1 - 2 * (k % 2))
            print(
                "Observer loop k = {} with Sun position "
                "at {} ".format(k + 1, pos_Sun)
            )

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
                    sdss_MV[select], sdss_MV[~select], sdss_class
                )
                sdss_bin[k * noSightings + i] = indiv_satellite_estimate(
                    pos2,
                    dis,
                    index,
                    sdss_surveyDirs[i, :],
                    sdss_cos_surveyAngle,
                    sdss_Rmax[select],
                    sdss_MV[select],
                    class_sats,
                    val_MV,
                )
                """ DES """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(des_MV.shape[0]) >= des_prob
                """ If any satellites are associated to LMC, update the
                classical satellite count """
                class_sats = update_classical_satellites(
                    des_MV[select], des_MV[~select], des_class
                )
                des_bin[k * noSightings + i] = indiv_satellite_estimate(
                    pos2,
                    dis,
                    index,
                    des_surveyDirs[i, :],
                    des_cos_surveyAngle,
                    des_Rmax[select],
                    des_MV[select],
                    class_sats,
                    val_MV,
                )
                """ SDSS + DES """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(tot_MV.shape[0]) >= tot_prob
                """ If any satellites are associated to LMC, update the
                classical satellite count """
                class_sats = update_classical_satellites(
                    tot_MV[select], tot_MV[~select], tot_class
                )
                tot_bin[k * noSightings + i] = combined_satellite_estimate(
                    pos2,
                    dis,
                    index,
                    sdss_surveyDirs[i, :],
                    sdss_cos_surveyAngle,
                    tot_Rmax_sdss[select],
                    des_surveyDirs[i, :],
                    des_cos_surveyAngle,
                    tot_Rmax_des[select],
                    tot_MV[select],
                    class_sats,
                    val_MV,
                )

        # Compile results
        SDSS_data_out[grp_i] = output_data_format(val_MV, sdss_bin, class_class)
        DES_data_out[grp_i] = output_data_format(val_MV, des_bin, class_class)
        SDSSDES_data_out[grp_i] = output_data_format(val_MV, tot_bin, class_class)

    # SDSS
    writeOutput(outputFile + "SDSS.hdf5", programOptionsDesc, SDSS_data_out, group_list)

    # DES
    writeOutput(outputFile + "DES.hdf5", programOptionsDesc, DES_data_out, group_list)

    # DES + SDSS
    writeOutput(
        outputFile + "SDSS+DES.hdf5", programOptionsDesc, SDSSDES_data_out, group_list
    )

    return None


# Run script
if __name__ == "__main__":
    main()
