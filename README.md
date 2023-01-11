# SPT-for-fMRI
MATLAB code for Shifted Partial Tracing for the analysis of fMRI

**Main programs:**

preprocess - given a motion-corrected nii file, perform the remainder of the preprocessing (cubic trend removal, cropping, standardisation)

pcaprojections3d - project data onto basis functions using pca

sptprojections3d - project data onto basis functions using spt

changepoints - given projected data, finds location of best fit AMOC and the corresponding test statistics

estdelta3d - computes the cost function for a range of values of bandwidth delta

analyse nii - given a motion-corrected nii, performs the remainder of the preprocessing, estimates bandwidth and performs change-point detection via PCA and SPT for a number of bandwidths

surveydata - given the 198 motion-corrected nii files and the txt file containing the list of subject names, performs change-point detection via PCA and SPT for a number of bandwidths

**Tools:**

einsum - https://uk.mathworks.com/matlabcentral/fileexchange/68995-einsum code by Yohai to perform a similar task to tensorprod

toeplitzouterprod, fastnormsquared, toeplitzsymbol - for the efficient computation of the cost function

**Data required:**

198 motion-corrected nii files

txt file containing list of subject names
