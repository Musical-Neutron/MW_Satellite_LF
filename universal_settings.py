#!/usr/bin/env python3

import os

# Cosmology
h = 0.73  # Hubble factor for the WMAP-1 cosmology used by Aquarius.

# File locations
all_sat_data = os.path.join('Input_Data', 'All_satellites.csv')
input_data = 'Input_Data'
output_data = 'Output_Data'

# Analysis settings
r_out = 300.
r_out_text = '300'
n_sightings = 1000
vmax_cut = 10.  # km/s

# Mass of the Aquarius haloes in units of 'Msun/h'.
aq_mass = {
    'Aq-A': 1.839e12 * h,
    'Aq-B': 0.819e12 * h,
    'Aq-C': 1.774e12 * h,
    'Aq-D': 1.774e12 * h,
    'Aq-E': 1.185e12 * h,
}

# A few values for the MW halo mass
MW_mass = {
    'small': [0.5e12 * h, '0.5e12'],
    'medium': [1.0e12 * h, '1.0e12'],
    'large': [2.0e12 * h, '2.0e12'],
}
