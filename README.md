# SPT-for-fMRI
MATLAB code for Shifted Partial Tracing for the analysis of fMRI
This was a summer research project completed by Alex Yan in 2022. Code is original, except for einsum.m which is from the MATLAB file exchange and carries its own license.

**Main programs:**

preprocess - given a motion-corrected nii file, perform the remainder of the preprocessing (cubic trend removal, cropping, standardisation)

pcaprojections3d - project data onto basis functions using pca

sptprojections3d - project data onto basis functions using spt

changepoints - given projected data, finds location of best fit epidemic change and the corresponding test statistics

estdelta3d - computes the cost function for a range of values of bandwidth delta

analyse_nii (script) - given a motion-corrected nii, performs the remainder of the preprocessing, estimates bandwidth and performs change-point detection via PCA and SPT for a number of bandwidths

**Tools:**

einsum - https://uk.mathworks.com/matlabcentral/fileexchange/68995-einsum code by Yohai (see einsum_license.txt) to perform a equivalent task to tensorprod

toeplitzouterprod, fastsquarednorm, fasttoeplitzsymbol - for the efficient computation of the cost function

**Not included:**

Motion-corrected nii files, unzipped
