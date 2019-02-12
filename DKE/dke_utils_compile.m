spm_path = 'Y:\helpern_data\Programs\spm8_r4667';  % path to spm 8, revision 4667
addpath(genpath(spm_path))

mcc  -v -m dke_preprocess_dicom.m ...
    -a spm_dicom_dict.mat

mcc  -v -m map_interpolate.m
