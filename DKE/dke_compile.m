spm_path = 'Y:\helpern_data\Programs\spm8_r4667';  % path to spm 8, revision 4667

addpath(genpath(spm_path))
addpath('C:\Users\Emilie\Documents\DKE\gradientVectors')

mcc  -m dke.m -a gradient_vectors_siemens6.dat -a gradient_vectors_siemens10.dat -a gradient_vectors_siemens12.dat -a gradient_vectors_siemens20.dat -a gradient_vectors_siemens30.dat -a gradient_vectors_siemens64.dat -a gradient_vectors_siemens256.dat

dke_utils_compile
