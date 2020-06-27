# SILVER - Set Increment with Limited Views Encoding Ratio
Pre-print available here: https://www.biorxiv.org/content/10.1101/2020.06.25.171017v1

## Setup
You can clone this project to get access to all code used in the SILVER project and re-run all experiments from the paper. If you want to add in pre-calculated SILVER optimizations, data and example results to save time, download [this](https://zenodo.org/record/3904738) dataset (~10GB) and add to the folder containing this repository and run setup_SILVER.m in Matlab. If you want to save space and setup time you can add only `SILVER_setup_basic.tar.gz` (about 500MB) from the dataset to your folder and run setup_SILVER.m. Then everything except example 7 (in-vivo) will run.

### Pre-requisites:
- [The NUFFT from the Michigan Image Reconstruction Toolbox](https://web.eecs.umich.edu/~fessler/code/)
- <b>Optional</b>: [MapVBVD](https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD)is needed to run  `generate_example_kdata` in `example7_invivo` if you want to use your own Siemens twix data.

## Use
### 2D
In 2D SILVER can be optimized with the Winkelmann method [1] or using electrostatic potential [2]. Simply call `ratio = SILVER_2D(window_sizes, method, save_filename)` with `window_sizes` being a list of integers representing how many spokes to combine in each window, `method` has to be a string that is either `Winkelmann` or `electrostatic_potential`, `save_filename` is optional and is the filename and path to where the result (`ratio` as well as other variables related to the optimistion (see comments in the function for details)) will be stored. `ratio` is the optimal step increment such that the angle to rotate each spoke is `ratio`*180 degrees. 

### 3D
Work in progress.

---

[1] Winkelmann S, Schaeffter T, Koehler T, Eggers H, Doessel O. An optimal radial profile order based on the Golden Ratio for time-resolved MRI. IEEE Trans Med Imaging. 2007;26(1):68-76. doi:10.1109/TMI.2006.885337

[2] Schauman SS, Okell TW, Chiew M. The Set Increment with Limited Views Encoding Ratio (SILVER) Method for Optimizing Radial Sampling of Dynamic MRI. BioRxiv 2020. https://doi.org/10.1101/2020.06.25.171017
