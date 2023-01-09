# SPT-for-fMRI
MATLAB code for Shifted Partial Tracing for the analysis of fMRI

Main programs:
analyse nii - given a motion-corrected nii, performs the remainder of the preprocessing, estimates bandwidth and performs change-point detection via (i) PCA (ii) SPT for a number of bandwidths
surveydata - given the 198 motion-corrected nii files and the txt file containing the list of subject names, performs change-point detection via SPT for a number of bandwidths

Subroutines:
nonsep3d - performs PCA to generate basis functions and executes changepoints on the projected data
spt3 - performs SPT to generate basis functions and executes changepoints on the projected data
changepoints - given projected data, finds location of best fit AMOC and the corresponding test statistics
einsum - https://uk.mathworks.com/matlabcentral/fileexchange/68995-einsum code by Yohai to perform a similar task to tensorprod
estdelta3d - computes the cost function for a range of values of bandwidth delta
toeplitzouterprod, fastnormsquared, toeplitzsymbol - for the efficient computation of the cost function

Data required:
198 motion-corrected nii files
txt file containing list of subject names
