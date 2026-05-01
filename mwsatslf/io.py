#!/usr/bin/env python3

import h5py
import numpy as np


def output_data_format(mv_bins, n_bins, class_class):
    """Prepares the data for writing.

    Args:
        mv_bins (arr): Magnitudes bins used in the estimation.
        n_bins (arr): Partially-complete array of estimated+classical
            luminosity functions.
        class_class (arr): Magnitudes of the classical satellites to
            prepend.

    Returns:
        arr (n_bins + 1): Luminosity functions. The first column
            corresponds to the absolute V-band magnitude values.
    """
    n_rows = mv_bins.shape[0] + class_class.shape[0]
    n_cols = 1 + n_bins.shape[0]
    out = np.empty((n_rows, n_cols), np.float32)

    # Insert the classical satellites
    class_class.sort()
    n_classical_sats = class_class.shape[0]
    # Magnitude (Mv) values corresponding to the classical satellites
    out[:n_classical_sats, 0] = class_class
    # Prepend cumulative classical satellite count
    out[:n_classical_sats, 1:] = (np.arange(n_classical_sats) + 1)[:, np.newaxis]

    # Append the estimated satellite count
    out[n_classical_sats:, 0] = mv_bins  # Magnitude bins
    out[n_classical_sats:, 1:] = n_bins.transpose()  # Estimated satellite count

    return out


def write_output(
    output_file,
    program_options_description,
    data_out,
    groups=False,
    n_sightings=1000,
    max_radius=300.0,
    min_vpeak=10.0,
):
    """Writes 'Lum_Func' data set to hdf5 file. This has the format
        returned by 'output_data_format'.

    Args:
        output_file (str): Full path to output file.
        program_options_description (str): Command line command used to produce
            this set of luminosity functions.
        data_out (tuple): The data to write to the file.
        groups (bool, optional): If list, list of groups to write into
            output file. Defaults to False.

    Returns:
        None
    """
    print(f"Writing the output data file '{output_file}' ...")
    with h5py.File(output_file, "w") as hf:
        if groups:
            for g_i, g_item in enumerate(groups):
                grp = hf.create_group(g_item)
                grp.create_dataset("Lum_Func", data=data_out[g_i], dtype="f4")
        else:
            hf.create_dataset("Lum_Func", data=data_out[0], dtype="f4")
        grp = hf.create_group("Metadata")
        grp.attrs["Program options"] = np.string_(program_options_description)
        n_options = program_options_description.split()
        if len(n_options) >= 7:
            grp.attrs["Number of sightings"] = n_options[6]
        else:
            grp.attrs["Number of sightings"] = n_sightings
        if len(n_options) >= 8:
            grp.attrs["Maximum radius"] = n_options[7]
        else:
            grp.attrs["Maximum radius"] = max_radius
        if len(n_options) >= 9:
            grp.attrs["Minimum V_peak"] = n_options[8]
        else:
            grp.attrs["Minimum V_peak"] = min_vpeak
    return None
