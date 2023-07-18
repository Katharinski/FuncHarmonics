### Information about the code contained in this repository

This repository is related to the publication "Functional harmonics reveal multi-dimensional basis functions 
underlying cortical organization", by K Glomb, ML Kringelbach, G Deco, P Hagmann, J Pearson, and S Atasoy, Cell Reports (2021) 

### Organization of the repository 

## code 

The pdf "link_code_analysis" contains a list of the analyses described in the paper and links them to the code. 

The main folder contains main scripts which perform the main parts of the analyses described in the paper. 

The folder "dependencies" contains functions called by these scripts as well as some scripts that need to be run once in order to produce certain supporting .mat files (like Matlab readable borders). The folder "thirdParty" contains the third party toolboxes (described below) for convenience. 


## data 

The folder "data/thisPaper" contains some of the data produced for the analyses described in the publication, namely: 

- the functional harmonics themselves (HCP_S900_CORR_manifold_knn300.mat) 
- the adjacency matrix with 300 nearest neighbors (adjacency_nn300.mat) 
- the cluster results necessary to reproduce some surface plots shown in Figure 2 (cluster_kmeans.mat) 

The folder "data/HCP_derived" contains files that we produced by converting HCP-provided files into different formats readable by Matlab, namely: 

- two label files derived from the HCP's RSN dense label files for the Yeo 7 RSN parcellation (RSN-networks.*.32k_fs_LR.label.gii, *=L or R) 
-border files in func.gii-format (see get_border_to_mat.m; borders_left.func.gii, boerder_right.func.gii, borders_som_left.func.gii borders_som_right.func.gii)

The folder "data/HCP" contains, for convenience, files downloaded from the HCP. They are described in the key resources table of the publication mentioned above. 

### Additional data needed to run this code 

FMRI data analyzed in the paper is not included in the repository due to its size. The key resources table of the publication mentioned above contains an exhaustive list of the data needed and where to obtain it. 

### 3rd party software needed to run this code 

As acknowledged in the paper mentioned above, we used the following 3rd party software packages included in this repository: 

- plot_mesh by Gabriel Peyr, https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/5355/versions/5/previews/toolbox_graph/plot_mesh.m/index.html (also includes getoptions.m and check_face_vertex.m)
- freezeColors by John Iversen, https://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors-unfreezecolors
- rgb.m by Ben Mitch, https://www.mathworks.com/matlabcentral/fileexchange/1805-rgb-m 
- subtightplot by Felipe G. Nievinski, https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
- MarkerTransparency toolbox by Free Software Foundation, https://www.mathworks.com/matlabcentral/fileexchange/65194-peterrochford-markertransparency

This is a static repository. An updated version of the code will be available on the corresponding author's github (github.com/Katharinski). Please send questions, bug reports and requests to katharina(dot)glomb(at)gmail(dot)com


