#!/usr/bin/env python3
import sys

import h5py
import numpy as np
from mwsatslf.geometry import uniform_points_on_sphere_surface
from mwsatslf.luminosity import (
    combined_satellite_estimate,
    indiv_satellite_estimate,
    update_classical_satellites,
)
from mwsatslf.io import output_data_format, write_output
from mwsatslf.cosmology import critical_density_200, compute_r200
from mwsatslf.survey import second_pointings, survey_Reff, survey_cone

survey_data = np.genfromtxt("Input_Data/Surveys.csv", names=True)
sat_data = np.genfromtxt("Input_Data/All_satellites.csv", names=True)
survey_area_dict = {
    "sdss": survey_data["SDSSIXW09"][0],
    "des": survey_data["DES"][0],
}  # Survey sky area
# Column label in satellite file that gives SDSS and DES satellites
dic_column_sats = {"sdss": "SDSSIX", "des": "DES"}

# a* and b* parameters to compute R_eff for a given magnitude.
dic_reff_params = {
    "sdss": {"a": survey_data["SDSSIXW09"][1], "b": survey_data["SDSSIXW09"][2]},
    "des": {"a": survey_data["DES"][1], "b": survey_data["DES"][2]},
}
n_sightings = 1000  # Default value for the number of sightings
mv_max = -8.8  # Brightest magnitude bin used for interpolation grid
n_bins_mv = 45
mv_min = mv_max * (n_bins_mv / (0.5 - n_bins_mv) + 1.0)
dmv = mv_max / (0.5 - n_bins_mv)
sun_distance = 8.0  # Distance of the sun from the Galactic Centre
distance_scaling_factor = 1.0e3  # Scale factor to obtain positions in 'kpc'
min_vpeak = 10.0  # Default value for the minimum V_peak used
max_radius = 400.0  # Default value for the maximum radius used
group_list = ["Original", "Original+Orphan", "Original+Orphan+Baryons"]


def main():
    program_name = sys.argv[0].rsplit("/")[-1]
    program_options_desc = program_name + " " + " ".join(sys.argv[1:])
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

    input_subhalo_file = sys.argv[1]
    host_mass = float(sys.argv[2])
    input_observation_file = sys.argv[3]
    output_host_mass = float(sys.argv[4])
    output_file = sys.argv[5]
    if len(sys.argv) >= 7:
        n_sightings = int(sys.argv[6])
    if len(sys.argv) >= 8:
        max_radius = float(sys.argv[7])
    if len(sys.argv) >= 9:
        min_vpeak = float(sys.argv[8])
    """Computes the value of '200 rho_critical' needed to compute R200
    given a halo mass M200. The numerical values correspond to the M200
    and R200 values of Aq.A1."""
    rho_crit_200 = critical_density_200(1.3432e12, 0.2458)

    # Obtain the fractional sky area/opening angle of each survey
    sdss_survey_area, sdss_cos_survey_angle = survey_cone(survey_area_dict["sdss"])
    des_survey_area, des_cos_survey_angle = survey_cone(survey_area_dict["des"])

    sdss_data_out = np.empty(len(group_list)).tolist()
    des_data_out = np.empty(len(group_list)).tolist()
    sdssdes_data_out = np.empty(len(group_list)).tolist()

    for grp_i, grp_item in enumerate(group_list):
        # Read in the subhaloes
        print("Reading subhalo data from file '{}' ...".format(input_subhalo_file))
        with h5py.File(input_subhalo_file, "r") as data:
            pos = np.asarray(data[grp_item + "/cop"])
            vpeak = np.asarray(data[grp_item + "/vpeak"])

        # Select only objects above the Vpeak cut
        s = vpeak >= min_vpeak
        s[0] = False  # Remove the host
        pos = pos[s] - pos[0]
        vpeak = vpeak[s]
        """ Compute a radial rescale factor. This is equivalent to changing
        the host mass to a new value. """
        # Input R200
        r200_original = compute_r200(host_mass, rho_crit200=rho_crit_200)
        # R200 for desired output halo mass
        r200_output = compute_r200(output_host_mass, rho_crit200=rho_crit_200)
        radius_scaling_factor = r200_output / r200_original
        print(
            "Rescaling subhalo positions by a factor of {:0.2f} "
            "corresponding to changing from the input mass, M200={:0.2e}, "
            "to the desired output mass, M200 = {:0.2e}.".format(
                radius_scaling_factor, host_mass, output_host_mass
            )
        )

        # Select only subhaloes within the desired distance from the host
        dis = (
            np.sqrt((pos * pos).sum(axis=1))
            * distance_scaling_factor
            * radius_scaling_factor
        )
        # Select everything with distances <= maximum radius + Solar distance
        s = dis <= (max_radius + sun_distance)
        pos = pos[s, :] * distance_scaling_factor * radius_scaling_factor
        vpeak = vpeak[s]

        # order the subhaloes according to their vMax values
        order = vpeak.argsort()[::-1]
        pos = pos[order, :]
        n_subhaloes = pos.shape[0]
        print(
            "There are {} subhaloes with Vmax >= {:0.1f} km/s and within a "
            "distance of {:0.0f} kpc from the central.".format(
                n_subhaloes, min_vpeak, max_radius + sun_distance
            )
        )
        """ Read the satellite data
        Reads magnitudes, distances, all provided satellite/survey data,
        and the probability of association with the LMC. """
        data = np.genfromtxt(input_observation_file, names=True)

        # SDSS satellites
        s = (data[dic_column_sats["sdss"]] == 2) * (data["Dkpc"] <= max_radius)
        mv_sdss = data["MV"][s]
        rad_sdss = data["Dkpc"][s]
        prob_sdss = data["LMC"][s]
        # Obtain classical satellites for SDSS
        s = (data[dic_column_sats["sdss"]] == 1) * (data["Dkpc"] <= max_radius)
        sdss_class = data["MV"][s]

        # DES satellites
        s = (data[dic_column_sats["des"]] == 2) * (data["Dkpc"] <= max_radius)
        mv_des = data["MV"][s]
        rad_des = data["Dkpc"][s]
        prob_des = data["LMC"][s]
        # Obtain classical satellites for DES
        s = (data[dic_column_sats["des"]] == 1) * (data["Dkpc"] <= max_radius)
        des_class = data["MV"][s]

        # Combine the two for the joint extrapolation
        mv_tot = np.append(mv_sdss, mv_des)
        rad_tot = np.append(rad_sdss, rad_des)
        prob_tot = np.append(prob_sdss, prob_des)
        # Obtain classical satellites for the combined SDSS and DES surveys
        s = (data[dic_column_sats["sdss"]] == 1) * (data["Dkpc"] <= max_radius) + (
            data[dic_column_sats["des"]] == 1
        ) * (data["Dkpc"] <= max_radius)
        tot_class = data["MV"][s]

        # Obtain 11 brightest classical satellites that are complete
        s = (data["Cl"] == 2) * (data["Dkpc"] <= max_radius)
        class_class = data["MV"][s]

        print(
            "There are {} satellites with magnitudes from {:0.1f} to "
            "{:0.1f}".format(mv_tot.shape[0], mv_tot.min(), mv_tot.max())
        )
        order = mv_sdss.argsort()
        mv_sdss = mv_sdss[order]
        prob_sdss = prob_sdss[order]

        order = mv_des.argsort()
        mv_des = mv_des[order]
        prob_des = prob_des[order]

        order = mv_tot.argsort()
        mv_tot = mv_tot[order]
        prob_tot = prob_tot[order]
        """ Generate a set of lines-of-sight (LOS) for the central direction
        of the survey, uniformly distributed over the surface of a sphere."""
        sdss_survey_directions = uniform_points_on_sphere_surface(n_sightings)
        des_survey_directions = second_pointings(sdss_survey_directions)

        # Compute the number of satellites using the 1st approach
        # Obtain magnitude bins
        bins_mv = np.linspace(mv_max, mv_min, n_bins_mv + 1, True)
        val_mv = 0.5 * (bins_mv[1:] + bins_mv[:-1])
        sdss_bin = np.zeros((6 * n_sightings, n_bins_mv), np.float32)
        des_bin = np.zeros((6 * n_sightings, n_bins_mv), np.float32)
        tot_bin = np.zeros((6 * n_sightings, n_bins_mv), np.float32)

        # Obtain bins for the individual satellites
        sdss_mv = np.append(bins_mv[0] - 0.1, np.append(mv_sdss, bins_mv[-1] + 0.1))
        sdss_rmax = survey_Reff(
            sdss_mv,
            a=dic_reff_params["sdss"]["a"],
            b=dic_reff_params["sdss"]["b"],
            survey_area=sdss_survey_area,
        )
        sdss_prob = np.append(0.0, np.append(prob_sdss, 0.0))

        des_mv = np.append(bins_mv[0], np.append(mv_des, bins_mv[-1] + 0.1))
        des_rmax = survey_Reff(
            des_mv,
            a=dic_reff_params["des"]["a"],
            b=dic_reff_params["des"]["b"],
            survey_area=des_survey_area,
        )
        des_prob = np.append(0.0, np.append(prob_des, 0.0))

        tot_mv = np.append(bins_mv[0], np.append(mv_tot, bins_mv[-1] + 0.1))
        tot_prob = np.append(0.0, np.append(prob_tot, 0.0))
        tot_rmax_sdss = survey_Reff(
            tot_mv,
            a=dic_reff_params["sdss"]["a"],
            b=dic_reff_params["sdss"]["b"],
            survey_area=sdss_survey_area,
        )
        tot_rmax_des = survey_Reff(
            tot_mv,
            a=dic_reff_params["des"]["a"],
            b=dic_reff_params["des"]["b"],
            survey_area=des_survey_area,
        )

        print(
            "Looping over 6 observer positions with {} sightings for "
            "each observer ...".format(n_sightings)
        )
        for k in range(6):  # Loop over observer positions
            # Obtain observer position
            sun_position = np.zeros(3, np.float32)
            """ For k=0, sets observer to x=+Sun_dis, k=1 -> x=-Sun_dis,
            k->2 y=+Sun_dis, and so on. """
            sun_position[int(k / 2)] = sun_distance * (1 - 2 * (k % 2))
            print(
                "Observer loop k = {} with Sun position "
                "at {} ".format(k + 1, sun_position)
            )

            # Translate subhalo positions relative to the observer position
            pos2 = pos + sun_position
            # Calculate distance of the subhaloes from the observer
            dis = np.sqrt((pos2 * pos2).sum(axis=1))
            select = dis <= max_radius

            # select the subhaloes within the required radius
            pos2 = pos2[select]
            dis = dis[select]
            index = np.arange(dis.shape[0])

            for i in range(n_sightings):  # Loop over the viewing directions
                # Randomises the subhalo list for each viewing direction.
                new_order = np.arange(len(pos2))
                np.random.shuffle(new_order)
                pos2 = pos2[new_order]
                dis = dis[new_order]
                """ SDSS """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(sdss_mv.shape[0]) >= sdss_prob
                """ If any satellites are associated to LMC, update the
            classical satellite count """
                class_sats = update_classical_satellites(
                    sdss_mv[select], sdss_mv[~select], sdss_class
                )
                sdss_bin[k * n_sightings + i] = indiv_satellite_estimate(
                    pos2,
                    dis,
                    index,
                    sdss_survey_directions[i, :],
                    sdss_cos_survey_angle,
                    sdss_rmax[select],
                    sdss_mv[select],
                    class_sats,
                    val_mv,
                )
                """ DES """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(des_mv.shape[0]) >= des_prob
                """ If any satellites are associated to LMC, update the
                classical satellite count """
                class_sats = update_classical_satellites(
                    des_mv[select], des_mv[~select], des_class
                )
                des_bin[k * n_sightings + i] = indiv_satellite_estimate(
                    pos2,
                    dis,
                    index,
                    des_survey_directions[i, :],
                    des_cos_survey_angle,
                    des_rmax[select],
                    des_mv[select],
                    class_sats,
                    val_mv,
                )
                """ SDSS + DES """
                # Select the subhaloes used for the extrapolation
                select = np.random.rand(tot_mv.shape[0]) >= tot_prob
                """ If any satellites are associated to LMC, update the
                classical satellite count """
                class_sats = update_classical_satellites(
                    tot_mv[select], tot_mv[~select], tot_class
                )
                tot_bin[k * n_sightings + i] = combined_satellite_estimate(
                    pos2,
                    dis,
                    index,
                    sdss_survey_directions[i, :],
                    sdss_cos_survey_angle,
                    tot_rmax_sdss[select],
                    des_survey_directions[i, :],
                    des_cos_survey_angle,
                    tot_rmax_des[select],
                    tot_mv[select],
                    class_sats,
                    val_mv,
                )

        # Compile results
        sdss_data_out[grp_i] = output_data_format(val_mv, sdss_bin, class_class)
        des_data_out[grp_i] = output_data_format(val_mv, des_bin, class_class)
        sdssdes_data_out[grp_i] = output_data_format(val_mv, tot_bin, class_class)

    # SDSS
    write_output(
        output_file + "SDSS.hdf5",
        program_options_desc,
        sdss_data_out,
        group_list,
        n_sightings,
        max_radius,
        min_vpeak,
    )

    # DES
    write_output(
        output_file + "DES.hdf5",
        program_options_desc,
        des_data_out,
        group_list,
        n_sightings,
        max_radius,
        min_vpeak,
    )

    # DES + SDSS
    write_output(
        output_file + "SDSS+DES.hdf5",
        program_options_desc,
        sdssdes_data_out,
        group_list,
        n_sightings,
        max_radius,
        min_vpeak,
    )

    return None


# Run script
if __name__ == "__main__":
    main()
