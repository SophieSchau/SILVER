# SILVER - Set Increment with Limited Views Encoding Ratio
New pre-print available here (currently under review): LINK HERE

Previous pre-print available here (rejected with invitation to resubmit when peer-reviewed): https://www.biorxiv.org/content/10.1101/2020.06.25.171017v1

## Abstract:
<b>Purpose</b>: To present and assess a method for choosing a fixed increment between spokes in radially sampled MRI that results in higher SNR than other radial sampling increments. 

<b>Theory and Methods</b>: Sampling uniformity contributes to image SNR when reconstructed using linear methods. Thus, for a radial trajectory, uniformly spaced sampling is ideal. However, uniform sampling lacks post-acquisition reconstruction flexibility, which is often needed in dynamic imaging. Golden ratio-based increments are often used for this purpose. The method presented here, Set Increment with Limited Views Encoding Ratio (SILVER), optimizes sampling uniformity when a limited number of temporal resolutions are required. With SILVER, an optimization algorithm finds the angular increment that provides the highest uniformity for a pre-defined set of reconstruction window sizes. SILVER was tested over multiple sets and assessed in terms of uniformity, noise amplification, and tSNR both in simulations and in acquisitions of a phantom and healthy volunteers. 

<b>Results</b>: The proposed algorithm produced trajectories that for the optimized window sizes had higher uniformity and lower image noise than golden ratio sampling both in a simulated single-coil system and in a multi-coil system, assessed using simulation, phantom, and in vivo experiments. The noise in SILVER optimized trajectories was comparable to uniformly distributed spokes whilst retaining flexibility in reconstruction at multiple temporal resolutions. In a resting state fMRI experiment, tSNR increases at different spatial/temporal resolutions were in the range 21-72%.

<b>Conclusion</b>: SILVER is a simple addition to any sequence currently using golden ratio sampling and significantly increases sampling efficiency and tSNR.


## Setup
You can clone this project to get access to all code used in the SILVER project and re-run all experiments. 


The figures from the paper can be re-generated using the script `paper2021.m` in the `experiments\` folder (for transparency, there is also a script for the rejected paper, `biorxiv2021.m`, but we do not recommend following it due to issues pointed out during peer review).

Data used in the experiments can be downloaded from LINK HERE (noise simulations), LINK HERE (k-space data for phantom and two volunteers), LINK HERE (reconstructions of phantom and volunteer data). The noise simulations were done using the function `SNR_predict_and_measure` The data was acquired on a 3T Siemens Prisma and the raw twix files were been pre-processed using the `read_in_and_save_TURBINE` function to generate the anonymised k-space data. Pre-calculated SILVER optimizations can also be downloaded from LINK HERE.

The data can be downloaded and put in appropriate directories by running the `download_data` script.

### Pre-requisites:
Must exist on the MATLAB search path:
- [The NUFFT from the Michigan Image Reconstruction Toolbox](https://web.eecs.umich.edu/~fessler/code/)

## Use
### 2D
In 2D SILVER can be optimized with the Winkelmann method [1] or using electrostatic potential [2] (recommended). Simply call `ratio = SILVER_2D(window_sizes, method, save_filename)` with `window_sizes` being a list of integers representing how many spokes to combine in each window, `method` has to be a string that is either `Winkelmann` or `electrostatic_potential`, `save_filename` is optional and is the filename and path to where the result (`ratio` as well as other variables related to the optimistion (see comments in the function for details)) will be stored. `ratio` is the optimal step increment such that the angle to rotate each spoke is `ratio`*180 degrees. 

### 3D
Work in progress.

---

[1] Winkelmann S, Schaeffter T, Koehler T, Eggers H, Doessel O. An optimal radial profile order based on the Golden Ratio for time-resolved MRI. IEEE Trans Med Imaging. 2007;26(1):68-76. doi:10.1109/TMI.2006.885337

[2] Schauman SS, Okell TW, Chiew M. The Set Increment with Limited Views Encoding Ratio (SILVER) Method for Optimizing Radial Sampling of Dynamic MRI. BioRxiv 2020. https://doi.org/10.1101/2020.06.25.171017
