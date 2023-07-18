%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for computing Euclidean distances between pairs of vertices on the
% functional harmonic and plotting it in an ordering that allows
% qualitative comparison with Yeo RSNs 
% 
% Copyright (c) 2004 Selen Atasoy
% Author: Selen Atasoy (selenatasoy@gmail.com)
% Date: 31/10/2018
%
% modified by Katharina Glomb, 2016
% katharina.glomb@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars

% load the manifold
load('HCP_S900_CORR_manifold_knn300.mat') % posCORR.M contains the functional harmonics

% load the surfaces
file_path = [pwd,'/data'];
subject_name = 'S900';
surf_type = 'inflated_MSMAll';
try
    [vertices, faces] = connRSMreadGII(file_path, subject_name, surf_type); % get surfaces
catch 
    error('Error reading surface. Perhaps gifti() function is not working. Try using gifti toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
end

% load the indices of the medial wall 
load('Ind_S900.mat')

save_now = false;

%% order the Glasser colormap
[color_ordering,glasser2016cmap_no_unknown] = get_color_order_HCP;
glasser2016cmap_final = [glasser2016cmap_no_unknown;glasser2016cmap_no_unknown];
nparc = length(glasser2016cmap_final);

%% load the surface labels
if ~exist('HCP_plot_labels_w_somatotopy.mat','file')
    error('File with surface labels of HCP parcellation including somatototopic areas does not exist. Get it by following instructions in create_som_subareas.m')
else
    load('HCP_plot_labels_w_somatotopy','surface_labels')
end
%% create the manifold
M = zeros(length(vertices.all), size(posCORR.M,2));
M(~(indices), :) = posCORR.M;

%% create the colormap
% the Glasser/HCP parcellation coloring follows an intricate scheme where
% blue-ness corresponds to the degree to which an area is related to visual
% functions, and equivalently for red:auditory, green:somatosensory/motor,
% white:task-positive, black:task-negative 
% this function returns parcels such that they are ordered
% auditory-somatosensory/motor-visual-taskNegative-taskPositive
% will be used for ordering within each RSN 
[inds_ordered] = [color_ordering;flipud(color_ordering+180)]; 

%% compute distances 
use_dims = 8:10; % fct harmonics used in paper: 8, 9, 10 (SI Figure...) 
counter = 0;

for d=use_dims
    counter = counter+1;
    
    fprintf(['Processing Dimension ', num2str(d), '...\n']);
    
    D = zeros(nparc,nparc);
    
    for m=1:nparc
        p1 = (surface_labels==m);
        D(m,m) = mean(pdist(M(p1,d+1), 'euclidean'));
        for n=m+1:nparc
            p2 = (surface_labels==n);
            D(m,n) = mean(mean(pdist2(M(p1,d+1), M(p2,d+1), 'euclidean')));
        end
    end
    D = D+D';
    M_norm = 1-(D/max(D(:)));
        
    % visualize RSN memberships 
    % Yeo colors 
    Visual = [120, 18, 134];
    Somatomotor = [70, 130, 180];
    Dorsal = [0, 118, 14];
    Ventral = [196, 58, 250];
    Limbic = [220, 248, 164];
    Frontoparietal = [230, 148, 34];
    Default = [205, 62, 78];
    Yeo_colors = [Visual; Somatomotor; Dorsal; Ventral; Limbic; Frontoparietal; Default]./255;
    cmap_clust_IDs = [Yeo_colors];%;rgb('black')];
    if ~exist('parcels_for_plotting_som.mat','file')
        error('File RSN memberships for HCP parcels is missing. Create it using main_connRSM_adjustYeoRSNs2HCP.m.')
    else
        load('parcels_for_plotting_som','parcel_IDs')
    end
    RSN_labels = zeros(nparc,1);
    for r=1:length(RSN_names{1})
        RSN_labels(parcel_IDs{1}{r}) = r;
    end
    [clust_labels,sort_ID] = sort(RSN_labels);
    
    % for clustering, also sort within the clusters according to color
    clust_labs_unique = unique(clust_labels);
    for c=1:length(clust_labs_unique)
        colorlabs_this_clust = sort_ID(clust_labels==clust_labs_unique(c));
        % find the subset of indices that correspond to this cluster
        target_this_clust = inds_ordered(ismember(inds_ordered,colorlabs_this_clust));
        % find each target index in the cluster indices to obtain ordering
        order_this_clust = zeros(length(target_this_clust),1);
        for e=1:length(target_this_clust)
            order_this_clust(e) = find(colorlabs_this_clust==target_this_clust(e));
        end
        sort_ID(clust_labels==clust_labs_unique(c)) = colorlabs_this_clust(order_this_clust);
    end

    %% actual plotting 
    % plot distance matrix itself
    twocolwidth = 18.3; % width of two col figures in cm
    h = figure('Units','centimeters','Position',[11.4375    5.8958    6.1648    4.4715],'PaperPositionMode','auto');
    subplot2 = @(m,n,p) subtightplot (m, n, p, [0 0], [0.01 0.01], [0.2 0.1]);
      
    % vectors that aide interpretation can be shown next to the matrix
    % HCP colors are always shown, plus one column/row of free space
    matrix_size = 17; 
    n_size = matrix_size+2; 
    % show RSN memberships in an extra row/column
    n_size = n_size+1;
    
    nvec = n_size-matrix_size;
    vec_pos_all = [];
    for v=nvec:-1:1
        vec_pos_left = v:n_size:(matrix_size*n_size);
        vec_pos_bottom = n_size*(n_size-v)+nvec+1:(n_size*(n_size-v+1));            
        if v==nvec-1 % vectors with HCP colors   
            % left
            subplot2(n_size,n_size,vec_pos_left);
            imagesc(sort_ID)
            axis off
            colormap(glasser2016cmap_final)
            freezeColors();
            % bottom
            subplot2(n_size,n_size,vec_pos_bottom);
            imagesc((sort_ID)')
            axis off
            colormap(glasser2016cmap_final)
            freezeColors();
        elseif v==nvec-2 
            % left 
            subplot2(n_size,n_size,vec_pos_left);
            imagesc(clust_labels)
            colormap(cmap_clust_IDs)
            axis off
            freezeColors();
            % bottom 
            subplot2(n_size,n_size,vec_pos_bottom);
            imagesc(clust_labels')
            axis off 
            colormap(cmap_clust_IDs)
            freezeColors();
        end
        vec_pos_all = [vec_pos_all,vec_pos_left,vec_pos_bottom];
    end
    
    % matrix itself
    matrix = setdiff(1:(n_size*n_size),vec_pos_all);
    subplot2(n_size,n_size,matrix(1:end-20));
    %axis square;
    imagesc(M_norm(sort_ID,sort_ID))
    colormap('gray')
    axis off;
    freezeColors();
    axis square
    set(gcf,'Color','white')   
    if save_now
        export_fig(['dist_mat_harmonic',  num2str(d), '.pdf'],'-r300');
    end
   
end