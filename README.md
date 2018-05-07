# Guidance for the MW satellite LF extrapolation scripts
**Last reviewed:** v1.0.0

**DOI:** 10.5281/zenodo.1205621

A code to calculate the luminosity function of Milky Way (MW) satellite galaxies
using Approximate Bayesian Computation. The analysis uses as input only the
observed satellite galaxies in the SDSS and DES surveys, and five DM-only
simulations of MW-like haloes from the Aquarius simulation suite (Aq-A – Aq-E).
Input data of the correct format is supplied for all known confirmed and
candidate satellite galaxies as of February 2018.

## 1.0 File and directory structure
This README file is located in the main repository for this analysis, from here
on referred to as the `home` directory. All required scripts can be found here.
Within this are the following directories:
```bash
.
├── Input_Data
│   ├── Aquarius
├── Output_Data

```

## 2.0 Simulation input data
The input data for each simulation are provided in the appropriate directory
under `home/Input_Data`, in the format:
   `<sim_name>_processed_trees.hdf5`.
Each group in this file corresponds to a different subhalo population, either:
* 'Original', the original population from the simulation
* 'Original+Orphan', with orphan galaxies identified and added as specified in
    [Newton+(2018)](https://doi.org/10.1093/mnras/sty1085)
* 'Original+Orphan+Baryons', the above but with baryonic correction.

For the [compute_satellites.py](/compute_satellites.py) script to work correctly
each group must include the following datasets, as a minimum:
* FOFCentre: 1 if the subhalo is the most massive subfind group in the FOF group.
* cop: The simulation x,y,z coordinates of the subhalo, in units of `Mpc`.
* vpeak: The peak maximum circular velocity of the subhalo attained at any time
         during its history, in units of `km/s`.

## 3.0 Satellite data file structure
The data file for all known MW satellite galaxies (including those which have
not yet been spectroscopically confirmed) is provided in the `home/Input_Data`
directory in `All_satellites.csv`.

In order, the columns of this file are as follows:
```
    'MV'       : Absolute V-band magnitude of the satelite   
    'Dkpc'     : Heliocentric distance of the satellite in kpc   
    'Cl'       : Flags indicating satellites to treat classically   
    'AT'       : Satellites in the VLT ATLAS survey   
    'DES'      : Satellites in the DES survey   
    'PAN'      : Satellites in the PanStaRRs survey   
    'SDSSV'    : Satellites in SDSS DR5   
    'SDSSIX'   : Satellites in SDSS DR9   
    'SDSSXII'  : Satellites in SDSS DR12   
    'SMASH'    : Satellites in the SMASH survey   
    'HSC'      : Satellites in the HyperSuprime Cam survey   
    'MagLiteS' : Satellites in the MagLiteS survey   
    'LMC'      : Probability that satellite is associated with LMC as   
                 determined by Jethwa et al. (2016) (http://bit.ly/2fymyJZ).   
    Last       : Names of the satellites
```

### 3.1 Survey Flags
    '0' excludes a satellite from analysis with a particular survey.   
    '1' in 'Cl' and another column: satellite will be treated classically when
        analysis is carried out with the specified survey.   
    '2' all other flags ignored. Satellite will be treated appropriately
        according to the survey it is in.   
        Note: '2' in the 'Cl' column indicates the 'original classicals'. These
        will appear in all analyses with any survey. This allows 'non-classical'
        satellites being treated classically to be distinguished from the 'true'
        ones.

The purpose of this structure is to allow, for example, some satellites on the
edge of SDSS DR5 to be treated classically, while in the expanded SDSS DR9
dataset they can be included.

## 4.0 Survey data file structure
The data files for all surveys that can be used in the analysis are provided in
the `home/Input_Data` directory in `Surveys.csv`. Each column corresponds to a
particular survey and response function. By default, this code ships with three:
```
    SDSSIX-K09 : SDSS DR9 using the Koposov+(2009) response function (http://bit.ly/1XaP5Qh)
    SDSSIX-W09 : SDSS DR9 using the Walsh+(2009) response function (http://bit.ly/2xAgSUy)
    DES        : DES using the Jethwa+(2016) response function (http://bit.ly/2fymyJZ)
```

Each column contains three rows:
```
    'Area' : Survey area in square degrees
    'a'    : Fitting parameter for the particular search algorithm, for use in
             Equation (3) of Newton+(2018) (https://doi.org/10.1093/mnras/sty1085).
    'b'    : Fitting parameter for the particular search algorithm, for use in
             Equation (3) of Newton+(2018) (https://doi.org/10.1093/mnras/sty1085).
```

## 5.0 Ouput data file structure
The `home/Output_Data` directory ships empty but can be filled with outputs from
[compute_satellites.py](/compute_satellites.py). By default the accompanying
[script](/script_compute_satellites.py) will write files to this location. These
will be in the format:
```
.
├── Original
│   ├── Lum_Func
├── Original+Orphan
│   ├── Lum_Func
├── Original+Orphan+Baryons
│   ├── Lum_Func
├── Metadata
```
* Lum_Func: An nxm array of luminosity functions, where n corresponds to the
number of magnitude bins (plus classical satellites), and m is the number of
luminosity functions produced per simulation, plus one column for Mv values.
* Metadata: Gives the command used to produce the data, the fiducial radius and
vpeak cut (in `kpc` and `km/s`, respectively), and the number of mock surveys
per observer position.

The groups correspond to those described in Section 2.0.

## 6.0 Citations
This code and accompanying input data are freely available. If using this code,
a derivative work or results thereof, please cite:
[Newton & Cautun (2017)](https://doi.org/10.5281/zenodo.1205621), and
[Newton+(2018)](https://doi.org/10.1093/mnras/sty1085)

If you have any questions or would like help in using the code, please email:
> oliver 'dot' j 'dot' newton 'at' durham.ac.uk
