function compileDKE
%   Compile DKE
%   Just run function to compile DKE into an executable for Mac
%   Ex:
%       >> compileDKE

%   Detect Platform
arch = computer('arch');
folderName = sprintf('DKE_Executable_%s',arch);

if ismac
    osName = 'Mac';
    mexFolder = ['.' filesep 'mac_mex'];
elseif ispc
    osName = 'Windows';
    mexFolder = ['.' filesep windows64_mex'];
elseif isunix
    osName = 'Unix';
    mexFolder = ['.' filesep 'linux_mex'];
end

fprintf('%s (%s) platform detected, building into "/%s"\n\n',...
        osName,arch,folderName);

%   Reset Relevant Folder
if isdir(folderName)
    rmdir(folderName);
    mkdir(folderName);
else
    mkdir(folderName);
end

%   Directories to include in compilation
gradientDir = ['.' filesep 'gradientVectors'];

%   Run compiler
fprintf('1. Compiling DKE (1/3)\n\n');
mcc('-v','-m','dke.m','-o','dke','-a',gradientDir,...
    '-a',mexFolder,'-d',folderName);

fprintf('2. Compiling dke_preprocess_dicom (2/3)\n\n');
mcc('-v','-m','dke_preprocess_dicom.m','-o','dke_preprocess_dicom',...
    '-a','spm_dicom_dict.mat','-d',folderName);

fprintf('3. Compiling map_interpolate (3/3)\n\n');
mcc('-v','-m','map_interpolate.m','-o','map_interpolation',...
    '-d',folderName);

fprintf('\n');
fprintf('completed\n\n')
end