#!/usr/bin/env python3
import sys
from subprocess import Popen
"""
Use it to execute a set of commands given as a list
Can be used with a batch system.

The script can be called without any arguments
or with two arguments, i.e.  name  N1  N2
where N1 and N2 are integers (starting with 0).
It will then execute only the commands in the
list starting from N1 up to N2.

e.g. ./script_lum_func_calculation.py 1 2
     will execute the 2nd and 3rd commands in
     the commands list.
"""

value = 1.0e12
value_text = '1.0e12'
h = 0.73  # Hubble factor for the WMAP-1 cosmology used by Aquarius.
Rout = 300.
Rout_text = '300'
# Mass of the Aquarius haloes in units of 'Msun/h'.
mass = {
    'Aq-A': 1.839e12 * h,
    'Aq-B': 0.819e12 * h,
    'Aq-C': 1.774e12 * h,
    'Aq-D': 1.774e12 * h,
    'Aq-E': 1.185e12 * h,
}

# A few values for the MW halo mass
MW_mass = {
    'small': 0.5e12 * h,
    'medium': 1.0e12 * h,
    'large': 2.0e12 * h,
}

commands = [
    # SMALL
    "./compute_satellites.py   Input_Data/Aquarius/Aq-A2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.A2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-A'], MW_mass['small'], '0.5e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-B2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.B2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-B'], MW_mass['small'], '0.5e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-C2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.C2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-C'], MW_mass['small'], '0.5e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-D2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.D2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-D'], MW_mass['small'], '0.5e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-E2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.E2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-E'], MW_mass['small'], '0.5e12', Rout_text, Rout),

    # MEDIUM
    "./compute_satellites.py   Input_Data/Aquarius/Aq-A2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.A2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-A'], MW_mass['medium'], '1.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-B2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.B2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-B'], MW_mass['medium'], '1.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-C2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.C2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-C'], MW_mass['medium'], '1.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-D2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.D2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-D'], MW_mass['medium'], '1.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-E2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.E2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-E'], MW_mass['medium'], '1.0e12', Rout_text, Rout),

    # LARGE
    "./compute_satellites.py   Input_Data/Aquarius/Aq-A2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.A2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-A'], MW_mass['large'], '2.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-B2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.B2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-B'], MW_mass['large'], '2.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-C2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.C2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-C'], MW_mass['large'], '2.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-D2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.D2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-D'], MW_mass['large'], '2.0e12', Rout_text, Rout),
    "./compute_satellites.py   Input_Data/Aquarius/Aq-E2_processed_trees.hdf5  {0}   Input_Data/All_satellites.csv  {1}  Output_Data/Aq.E2_M_MW_{2}_r_{3}_ 1000 {4} 10."
    .format(mass['Aq-E'], MW_mass['large'], '2.0e12', Rout_text, Rout),
]

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
