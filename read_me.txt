Thank you for your interest in FBWM! 

The main code can be found under Scripts/FBWM.m and requires the following inputs: 

- Path_data: path to a 4D nifti file with preprocessed(*) diffusion data. 
- Path_gradient: path to a txt file with the matching diffusion gradients (Nx3).
- Path_bval: path to a text file with the matching b-values (Nx1).  
- Path_brain_mask: path to a binary brain mask. 
- Path_output: path to where the output will be written.
- path_dke_parameters: path to template DKE parameters file that was provided with code (don't make any changes to this file). 
- degree: degree of spherical harmonics used to estimate the fODF (max degree = 8) 

Before running add the DKE folder to your Matlab path. Requires SPM to read and write niftis. 

An example dataset with corresponding outputs is located in the Example folder. 

When using this code please cite the following papers: 
McKinnon, E. T., Helpern, J. A., & Jensen, J. H. (2018). Modeling white matter microstructure with fiber ball imaging. NeuroImage, 176, 11-21.
Jensen, J. H., Glenn, G. R., & Helpern, J. A. (2016). Fiber ball imaging. Neuroimage, 124, 824-833.


* Before running FBWM, we always run at a minimum the following preprocessing steps:
- MP-PCA denoising ( https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html ) 
- Gibbs Artefact removal (https://mrtrix.readthedocs.io/en/latest/reference/commands/mrdegibbs.html )  
- Rician bias correction 
- Motion correction  

Please refer to the following reference for more details on dMRI preprocessing: Ades-Aron, Benjamin, et al. "Evaluation of the accuracy and precision of the diffusion parameter EStImation with Gibbs and NoisE removal pipeline." NeuroImage (2018).

Feel free to contact us with any questions or remarks (mckinnon@musc.edu). 
