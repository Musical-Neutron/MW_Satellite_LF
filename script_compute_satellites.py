#!/usr/bin/env python3
import os
import sys
from subprocess import Popen
from universal_settings import (r_out, r_out_text, aq_mass, MW_mass,
                                all_sat_data, input_data, output_data)
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

# value = 1.0e12
# value_text = '1.0e12'
# h = 0.73  # Hubble factor for the WMAP-1 cosmology used by Aquarius.
# r_out = 300.
# r_out_text = '300'
# # Mass of the Aquarius haloes in units of 'Msun/h'.
# mass = {
#     'Aq-A': 1.839e12 * h,
#     'Aq-B': 0.819e12 * h,
#     'Aq-C': 1.774e12 * h,
#     'Aq-D': 1.774e12 * h,
#     'Aq-E': 1.185e12 * h,
# }

# # A few values for the MW halo mass
# MW_mass = {
#     'small': 0.5e12 * h,
#     'medium': 1.0e12 * h,
#     'large': 2.0e12 * h,
# }

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
    'Aq.A2_M_MW_{0}_r_{1}_', 'Aq.B2_M_MW_{0}_r_{1}_', 'Aq.C2_M_MW_{0}_r_{1}_',
    'Aq.D2_M_MW_{0}_r_{1}_', 'Aq.E2_M_MW_{0}_r_{1}_'
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
        command = ('./compute_satellites.py ' + '{} '.format(tree) +
                   '{} '.format(aq_mass[aq_identifier]) + all_sat_data +
                   ' {}'.format(MW_mass[mw_size][0]) + ' {}'.format(
                       output_file.format(MW_mass[mw_size][1], r_out_text) +
                       ' 1000' + ' {}'.format(r_out) + ' 10.'))
        commands.append(command)

# commands = [
#     # SMALL
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-A2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.A2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-A'], MW_mass['small'][0], MW_mass['small'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-B2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.B2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-B'], MW_mass['small'][0], MW_mass['small'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-C2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.C2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-C'], MW_mass['small'][0], MW_mass['small'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-D2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.D2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-D'], MW_mass['small'][0], MW_mass['small'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-E2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.E2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-E'], MW_mass['small'][0], MW_mass['small'][1],
#             r_out_text, r_out),

#     # MEDIUM
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-A2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.A2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-A'], MW_mass['medium'][0], MW_mass['medium'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-B2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.B2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-B'], MW_mass['medium'][0], MW_mass['medium'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-C2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.C2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-C'], MW_mass['medium'][0], MW_mass['medium'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-D2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.D2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-D'], MW_mass['medium'][0], MW_mass['medium'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-E2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.E2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-E'], MW_mass['medium'][0], MW_mass['medium'][1],
#             r_out_text, r_out),

#     # LARGE
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-A2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.A2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-A'], MW_mass['large'][0], MW_mass['large'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-B2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.B2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-B'], MW_mass['large'][0], MW_mass['large'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-C2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.C2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-C'], MW_mass['large'][0], MW_mass['large'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-D2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.D2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-D'], MW_mass['large'][0], MW_mass['large'][1],
#             r_out_text, r_out),
#     "./compute_satellites.py   Input_Data/Aquarius/Aq-E2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.E2_M_MW_{2}_r_{3}_ 1000 {4} 10."
#     .format(aq_mass['Aq-E'], MW_mass['large'][0], MW_mass['large'][1],
#             r_out_text, r_out),
# ]

#
# No need to modify below this
#

noParts = len(commands)
minPart, maxPart = 0, noParts
if len(sys.argv) >= 3:
    minPart, maxPart = int(sys.argv[1]), int(sys.argv[2]) + 1
    if minPart > noParts: sys.exit(1)
    if maxPart > noParts: maxPart = noParts

# run the code for each part
procs = []
for i in range(minPart, maxPart):  # loop over all parts
    toRun = commands[i]
    print("\nRUNNING:\n\{}\n".format(toRun))
    proc = Popen(toRun, shell=True)
    procs.append(proc)

for proc in procs:
    proc.wait()
