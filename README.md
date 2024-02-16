# Scattering Horizons

This repository provides Python scripts that simulate the DM and scattering time (in ms at 1 GHz) contributions from galaxies along FRB lines-of-sight (LOSs), including nearby galaxies in the Gravitational Wave Catalog (GWGC), more distant intervening galaxies (modeled statistically using a halo mass function), and FRB host galaxies. The Milky Way can be modeled separately using NE2001 or YMW16 and an electron density prescription for the Galactic halo. Please cite [Ocker et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...934...71O/abstract) for use of this repository. A complete description of the methods can be found in the paper.

## Dependencies and Installation

The provided scripts have the usual Python dependencies (astropy, numpy, scipy, matplotlib). Additionally, the colossus package is required and can be downloaded from pip (`pip install colossus`) or from source code (see [package website](https://pypi.org/project/colossus/)). 

All of the provided scripts should be downloaded into the same working directory, as several of the scripts provide functions that are imported for use in the other scripts.

## Code Description 

There are five scripts in the repository:

**tauDMdisk.py**

Evaluates the scattering time and DM for a single LOS through a galaxy ISM. Simple electron density prescriptions are provided for three galaxy types, dwarfs, ellipticals, and spirals. Functions defined here are imported to intervening_galaxies.py and host_galaxies.py to simulate many LOSs.

**intersection_probabilities.py**

Evaluates the number of galaxies along the LOS to a given source redshift, using the Tinker halo mass function. This script is imported into intervening_galaxies.py.

**gwgc_tau.py**

Predicts DM and scattering times from galaxies in GWGC. The script reads in the GWGC catalog ('gwgc_cat.fit'). 

*To use this code:* 
- Modify the value of the FRB source redshift, 'frb_zs' on Line 86.
- Modify the path to the output file on Line 197. The output is saved as a .npz file that can be read using the numpy load() routine.

**intervening_galaxies.py**

Simulates DMs and scattering times from intervening galaxies at moderate redshift (i.e., beyond the Local Volume captured in gwgc_tau.py). Galaxy types are defined by stellar mass, which is estimated based on the stellar-to-halo mass function from [Girelli et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020A%26A...634A.135G/abstract). Please cite Girelli et al. if you make use of this script.

*To use this code:* 
- Modify the value of the FRB source redshift on Line 202. 
- Modify the redshift sampling interval defined in the redshift array 'z_arr' on Line 88; this sampling interval determines the number of redshift bins used by intersection_probabilities.py to estimate the number of intervening galaxies. E.g., for an FRB redshift of 5, I used a sampling interval of 1 whereas for an FRB redshift of 1, I used a sampling interval of 0.2 (i.e., I went for 5 redshift bins, mainly for speed).
- You can also modify the number of FRBs to simulate by changing the shapes of 'tautot_arr' and 'dmtot_arr' on Lines 208-209 (the default setup is to simulate 361x181 uniformly distributed sightlines, which takes awhile). 
- Be sure to modify the output file path on Line 230 (last line of the script).

**host_galaxies.py**

Simulates host galaxy contributions, assuming FRBs are uniformly distributed within their hosts (i.e., no specific formation channel is assumed). In addition to DMs and scattering times, the code outputs the host galaxy type given to each FRB. Note that the host galaxy type is assigned probabilistically based on the the galaxy stellar mass function from [McLeod et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.4413M/abstract). Please cite McLeod et al. if you make use of this script.

*To use this code:*
- DM and scattering times are generated from Line 306 onwards for a few different source redshifts, by evaluating the function 'host_galaxies_eval(frb_zs)'. Modify the value of 'frb_zs' and or comment out unnecessary lines of code.
- Modify the output file path.
- Modify the number of FRBs to simulate by changing the shapes of 'tautot_arr', 'dmtot_arr', and 'galtype_arr' on Lines 214-216. The default is to simulate 361x181 FRBs.
- Optionally, you can change the maximum radius within which FRBs are produced. The default is to assume FRBs can occur anywhere within the galaxy disk (scaled to 19 kpc for a Milky Way-mass spiral galaxy); the relevant variable to change is 'rdisk' defined on Lines 270, 282, and 294 (in three separate code blocks for the three different galaxy types). 

## Some Caveats

Halo gas is assumed to extend out to twice the virial radius, which can substantially increase the DM relative to simply using the virial radius. The impact on scattering is less straightforward (depending on halo mass, the corresponding increase in path length can counteract the increase in DM; recall tau~DM^2/L). You can modify the halo cutoff under the function 'ne_halo()' in the above scripts (host_galaxies, intervening_galaxies, gwgc_tau): change 'necut = ne_mnfw*0.5*(1-tanh((r-2.*r200)/20.))' to 'necut = ne_mnfw*0.5*(1-tanh((r-r200)/20.))'.

Dwarf galaxies dominate the intervening galaxy distribution because they are so numerous, but this means that the DM distribution presented in Ocker et al. (2022) appears to peak at smaller DM values than in some other studies based on full cosmological simulations. If the DM contributions of intervening galaxies is being underestimated, then the scattering may also be underestimated. In principle, one could modify the halo mass function used to draw the intervening galaxy masses, via the intersection_probabilities.py script.
