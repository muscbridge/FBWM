thisDirectory = pwd;

[rv, errorMsg] = mkdir('win64');
if rv ~= 1
    fprintf(errorMsg)
    return
end

% [rv,msg] = copyfile('mfiles\*.m', 'win64');
% if rv ~= 1
%     fprintf(errorMsg)
%     return
% end

spmPath = 'D:\spm8'; 
addpath(genpath(spmPath))

mexPath = [thisDirectory, '\win64mex'];
addpath(genpath(mexPath))

% cd('win64')
cd('mfiles')

% gradientVector_path = '.\gradientVectors';
% addpath(genpath(gradientVector_path))

% mcc -v -m dke.m ...
%     -a gradient_vectors_siemens6.dat ...
%     -a gradient_vectors_siemens10.dat ...
%     -a gradient_vectors_siemens12.dat ...
%     -a gradient_vectors_siemens20.dat ...
%     -a gradient_vectors_siemens30.dat ...
%     -a gradient_vectors_siemens64.dat ...
%     -a gradient_vectors_siemens256.dat

% mcc -v -m dke.m -a ..\gradientVectors\*.dat
% 
% % dke_utils_compile
% mcc -m dke_preprocess_dicom.m -a spm_dicom_dict.mat
% mcc -m map_interpolate.m

outputPath = [thisDirectory, '\win64'];

fprintf('building DKE (win64, external)\n')

fprintf('compiling DKE\n')
mcc -d ..\win64 -o dke -m dke.m -a ..\gradientVectors\*.dat

fprintf('compiling dke_preprocess_dicom\n')
mcc -d ..\win64 -o dke_preprocess_dicom -m dke_preprocess_dicom.m -a spm_dicom_dict.mat

fprintf('compiling map_interpolate\n')
mcc -d ..\win64 -o map_interpolate -m map_interpolate.m

fprintf('completed\n\n')

cd(thisDirectory)


% cd('win64_internal')
% dke_compile
% cd(thisDirectory)

