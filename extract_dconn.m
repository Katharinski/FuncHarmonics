% script to extract vector of unique cortical pairs from HCP file, and 
% save in .mat-format
% Katharina Glomb 
% last updated: June 22, 2021 
% katharina.glomb@gmail.com 

%%%% important information %%%%
% --loading the HCP dense connectome requires extra memory
% the file contains a single variable called cdata, which is a 91282 x 91282
% single precision matrix 
% --for loading this kind of file (cifti format), HCP provides a toolbox
% see https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ
% Question 2 "How do I get cifti files into Matlab?" 
% We tested all options with the same results. 
% --apart from this, connectome workbench is necessary (see the same FAQ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars

% the paths assume Linux 
% applications are installed in "snap"
% code of this repository is in "repos" and the repo is named "FuncHarmonics"
% you placed the group FC in repos/FuncHarmonics/data/HCP

homedir = getenv('HOME');
repodatadir = fullfile(homedir,'repos','FuncHarmonics','data');

% connectome workbench
wb_dir = fullfile(homedir,'snap','workbench/bin_linux64/wb_command');
FCfile = fullfile(repodatadir,'HCP','HCP_S1200_812_rfMRI_MSMAll_groupPCA_d4500ROW_zcorr_recon2.dconn.nii');
ciftiData = ciftiopen(FCfile,wb_dir);

% get cortical indices; this contains the medial wall, but dconn DOES NOT
load(fullfile(repodatadir,'HCP_derived','Ind_S900'),'indices')

nvox = sum(~indices); % #non-medial wall voxels

cdata = ciftiData.cdata(1:nvox,1:nvox);

clear ciftiData

% turn cdata into a distance vector
cdata = 1-tanh(cdata); % matrix contains z-transformed correlations
cdata(logical(eye(size(cdata)))) = 0; % make sure diagonal is 0
cdata = squareform(cdata);
save('/mnt/data/FuncHarmonics/data/HCP_S900_pdist_vec','cdata','-v7.3')
