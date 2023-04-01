#!/usr/bin/env python3
import os
from subprocess import Popen

from universal_settings import (MW_mass, all_sat_data, aq_mass, input_data,
                                n_sightings, output_data, r_out, r_out_text,
                                vmax_cut)
"""
Use it to execute a set of commands given as a list
Can be used with a batch system.

The script can be called without any arguments
or with two arguments, i.e.  name  N1  N2
where N1 and N2 are integers (starting with 0).
It will then execute only the commands in the
list starting from N1 up to N2.

e.g. ./script_compute_satellites.py 1 2
     will execute the 2nd and 3rd commands in
     the commands list.
"""


def main():
    input_tree_files = [
        'Aq-A2_processed_trees.hdf5', 'Aq-B2_processed_trees.hdf5',
        'Aq-C2_processed_trees.hdf5', 'Aq-D2_processed_trees.hdf5',
        'Aq-E2_processed_trees.hdf5'
    ]
    input_trees = [
        os.path.join(input_data, 'Aquarius', tree_file)
        for tree_file in input_tree_files
    ]
    output_file_templates = [
        'Aq.A2_M_MW_{0}_r_{1}_', 'Aq.B2_M_MW_{0}_r_{1}_',
        'Aq.C2_M_MW_{0}_r_{1}_', 'Aq.D2_M_MW_{0}_r_{1}_',
        'Aq.E2_M_MW_{0}_r_{1}_'
    ]
    output_files = [
        os.path.join(output_data, file_template)
        for file_template in output_file_templates
    ]
    mw_size_list = list(MW_mass.keys())

    commands = []
    for tree, output_file in zip(input_trees, output_files):
        aq_identifier = tree.split('/')[-1].split('_')[0][:4]
        for mw_size in mw_size_list:
            command = (
                './compute_satellites.py ' + '{} '.format(tree) +
                '{} '.format(aq_mass[aq_identifier]) + all_sat_data +
                ' {}'.format(MW_mass[mw_size][0]) + ' {}'.format(
                    output_file.format(MW_mass[mw_size][1], r_out_text) +
                    ' {}'.format(int(n_sightings)) + ' {}'.format(r_out) +
                    ' {}'.format(vmax_cut)))
            commands.append(command)

    #
    # No need to modify below this
    #

    # run the code for each generated command
    procs = []
    for toRun in commands:  # loop over all parts
        print("\nRUNNING:\n\{}\n".format(toRun))
        proc = Popen(toRun, shell=True)
        procs.append(proc)

    for proc in procs:
        proc.wait()

    return None


if __name__ == "__main__":
    main()
