   %% "FBI"Äù: Fiber Ball Imaging 
    %       input:
    %           - path_data: path to 4D nifti 
    %           - path_gradient:   path to gradient directions in same order as data
    %           - path_bval:  path to bvalues in same order as data
    %           - degree: degree of spherical harmonics used for fODF 
    %           - path_brain_mask: path to binary mask 
    %           - path_output: path to output data
    %           - path_dke_parameters= path to template file DKEParameters.txt

    %   Authors: Russell Glenn, Hunter Moss and Emilie McKinnon (mckinnon@musc.edu), Updated August 31, 2018
    %   Copyright (c) 2018 medical university of south carolina (MUSC)
    %       
    %      Permission is hereby granted, free of charge, to Prof. Wheeler-Kingshott's lab
    %      to use the Software solely for non-commercial research, including the rights to use
    %      and modify the Software, subject to the following conditions:
    %       
    %        1. The above copyright notice and this permission notice shall be
    %      included by Recipient in all copies or substantial portions of the
    %      Software. 
    %       
    %        2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
    %      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF
    %      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    %      NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
    %      DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    %      OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
    %      USE OR OTHER DEALINGS IN THE SOFTWARE. 
    %       
    %        3. In no event shall MUSC be liable for direct, indirect, special,
    %      incidental or consequential damages in connection with the Software.
    %      Recipient will defend, indemnify and hold MUSC harmless from any claims or
    %      liability resulting from the use of the Software by recipient. 
    %       
    %        4. Neither anything contained herein nor the delivery of the Software to
    %      recipient shall be deemed to grant the Recipient any right or licenses
    %      under any patents or patent application owned by MUSC. 
    %       
    %        5. The Software may only be used for non-commercial research and may not
    %      be used for clinical care. 
    %       
    %        6. Any publication by Recipient of research involving the Software shall
    %      cite the references listed below.
    %
    %     REFERENCES
    %     Jensen, J. H., Glenn, G. R., & Helpern, J. A. (2016). Fiber ball imaging. Neuroimage, 124, 824-833.
 


%% input 

path_data='';
path_gradient='';
path_brain_mask='';
degree=6;
path_bval='';
path_output='';

%%
% read in bvalues and gradient table 
bt= textscan(fopen(fullfile(path_bval)),'%f');
gt= textscan(fopen(fullfile(path_gradient)),'%f%f%f');

% define bvalues used
bval=unique(bt{1});
bval_fbi=bval(end)
index_fbi= find(bt{1}==bval_fbi);

% number of gradient directions for each b-value 
ndir_fbi=length(index_fbi);

% check units of bvalues and put it in ms/um2 
if (bval_fbi/1000)<1
else
bval_fbi=bval_fbi/1000; 
end

%Separate b0 from DWIs
idx_b0 = find(bt{1}==0); 

% Read in data
S=spm_read_vols(spm_vol(fullfile(path_data)));
dim=size(S);

Sb6_reshape=zeros(1,prod(dim(1:3)));

S_b6=S(:,:,:,index_fbi);
for i=1:length(index_fbi)
    Sb6_reshape(i,:)=reshape(S_b6(:,:,:,i),[1,prod(dim(1:3))]);
end

S_b0=S(:,:,:,idx_b0);
S_b0_mean=mean(S_b0,4); 
s0=reshape(S_b0_mean,[1, prod(dim(1:3))]);

% brain mask 
hdr=spm_vol(path_brain_mask);
brainmask=spm_read_vols(hdr);
brain_mask_reshape=reshape(brainmask,[1,prod(dim(1:3))]);

% seperate gradient tables 
GT_fbi = [gt{1}(index_fbi) gt{2}(index_fbi) gt{3}(index_fbi)];  
GT_fbi=GT_fbi./sqrt(GT_fbi(:,1).^2+GT_fbi(:,2).^2+GT_fbi(:,3).^2);

GT={GT_fbi};

%GET SH BASIS FUNCTIONS
B_fbi=getSH(degree,[ atan2(GT_fbi(:,2),GT_fbi(:,1)) acos(GT_fbi(:,3))],'complex');
Harm_id={1,5:9,17:25,37:49,65:81}; % max degree = 8 
B_fbi=B_fbi(:,cell2mat(Harm_id(1:degree/2+1)));
B={B_fbi};


%% Core calculations (can take about 10 min) 
% initialisation 
tic 
intermediate_A_img=zeros(3,3,prod(prod(dim(1:3))));
intermediate_DT_img=zeros(3,3,prod(prod(dim(1:3))));
zeta_img=zeros(1,prod(dim(1:3)));
faa_f=zeros(1,prod(dim(1:3)));

% calculate legendre functions and g2l 
D0=3;
idx_Y = 0;
    for l = 0:degree/2
        P2l0(idx_Y+1:idx_Y+(4*l+1), :) =  ((-1)^(l) * factorial(2*l))/(4^(l)*factorial(l)^2)*ones(4*l+1,1); %equation 12 FBI paper
        G2l(idx_Y+1:idx_Y+(4*l+1), :) = (factorial(l)*(bval_fbi*D0)^((l+0.5))/gamma(2*l+3/2)*hypergeom((l+0.5),2*l+3/2,-bval_fbi*D0))*ones(4*l+1,1); %equation18 FBI paper
        idx_Y = idx_Y + 4*l+1;
    end
    
% voxel-wise estimation of other parameters 
parfor vox=1:prod(dim(1:3))

%voxelwise signal values
S0_vox=s0(vox);
Sfbi_vox=Sb6_reshape(:,vox);
S={Sfbi_vox};
    if brain_mask_reshape(vox)==1

     
A2l=(B_fbi'*B_fbi)^-1*B_fbi'*Sfbi_vox/s0(vox); %equation 4 FBWM paper 
C2l=A2l*G2l(1).*(sqrt(4*pi)*P2l0*A2l(1).*G2l).^-1; %equation 8 FBWM paper 

C2l = bsxfun(@rdivide,C2l,C2l(1)); %normalization
C2l= bsxfun(@times,C2l,1/(sqrt(4*pi))); %normalization
zeta=A2l(1)*sqrt(bval_fbi)./pi;

scale=repmat((((C2l(1)*sqrt(30))).^-1),[1 6]); % scale factor A 

A=scale.*[sqrt(30)/3*C2l(1)-sqrt(6)/3*C2l(4)+C2l(6)+C2l(2),...
            sqrt(30)/3*C2l(1)-sqrt(6)/3*C2l(4)-C2l(6)-C2l(2),...
            sqrt(30)/3*C2l(1)+2*sqrt(6)/3*C2l(4),...
            1i*C2l(6)-1i*C2l(2),...
            -C2l(5)+C2l(3),...
            -1i*C2l(5)-1i*C2l(3)]; %equation 12 FBWM paper 

intermediate_A=[A(1), A(4), A(5);
                A(4), A(2), A(6);
                A(5), A(6), A(3)];

faa_f(vox)= sqrt(3*sum(abs(C2l(2:6)').^2)./(5*abs(C2l(1)').^2 + 2*sum(abs(C2l(2:6)').^2))); % equation 13 FBWM 
zeta_img(vox)=zeta;
intermediate_A_img(:,:,vox)=intermediate_A;
    end
    
end
toc

%% write outputs 
mkdir([path_output '/FBI/']);
hdr.dt=[16 0];
zeta_reshape=reshape(zeta_img,dim(1:3));hdr.fname=[path_output '/FBI/zeta.nii'];spm_write_vol(hdr,zeta_reshape);
faa_f_reshape=reshape(faa_f,dim(1:3));hdr.fname=[path_output '/FBI/faa_f.nii'];spm_write_vol(hdr,faa_f_reshape);

save([path_output '/FBI/A.mat'],'intermediate_A_img');
save([path_output '/FBI/intermediate_A_img.mat'],'intermediate_A_img');
