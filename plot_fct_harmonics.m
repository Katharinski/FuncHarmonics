% produce plots of functional harmonics on the cortical surface (figure 2
% of the paper), adding borders 
% code by Selen Atasoy & Katharina Glomb 
% last modified: June 26, 2021
% katharina.glomb@gmail.com 

clearvars

homedir = getenv('HOME');
repodatadir = fullfile(homedir,'repos','FuncHarmonics','data');

load('Ind_S900','indices') % medial wall indices; see extract_dconn.m
nvox = sum(~indices); % #non-medial wall voxels

% load cortical surfaces 
path_to_surfaces = fullfile(repodatadir,'HCP');
subject_name = 'S900';
surf_type = 'inflated_MSMAll';
% surf_type = 'flat';
try
    [vertices, faces] = connRSMreadGII(path_to_surfaces, subject_name, surf_type); % get surfaces
catch
    error('Error reading surface. Perhaps the path to the surface file is incorrect, or gifti() function is not working. Try using gifti toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
end

%% some parameters
border_color = 'white'; % 'none'

% adjust number of subplots depending on which surface is used and whether
% you want a colorbar or not (inflated: 6 views, flat: 2 views) 
plot_colorbar = true;
if strcmp(surf_type,'inflated_MSMAll') && plot_colorbar
    nsub = 7;
elseif strcmp(surf_type,'inflated_MSMAll') && ~plot_colorbar
    nsub = 6;
elseif strcmp(surf_type,'flat') && plot_colorbar
    nsub = 3;
elseif strcmp(surf_type,'flat') && ~plot_colorbar
    nsub = 2;
end
%% plot
% options for subtightplot
gapleft = 0.035;
gapright = 0.035;
margin_height = [0.03 0.03]; % lower, upper
margin_height_flat = [0.15 0.15]; % lower, upper
margin_width = [0.04 0.06]; % left, right

load('HCP_S900_CORR_manifold_knn300.mat','posCORR') % posCORR.M contains fct harmonics
% note: fct harmonics are multiplied *1000 to make colorbars manageable;
% the label "x 10^-3" is added manually to the plots as Matlab cannot
% handle this 
for mf=2 % mf refers to the column of posCORR.M; first non-constant fct harmonic: mf=2
    f = figure;
    % for some harmonics, we show k-means clusters; manually say which clusters to plot
    choose_clusters = zeros(12,2);
    choose_clusters([6,7,9,11],1) = 1; 
    choose_clusters([6,7,9,11],2) = 2;    
    choose_clusters(11,1) = 3;
    
    % add medial wall vertices (which are not present in posCORR.M) 
    plot_data_mf = zeros(length(indices),1);
    plot_data_mf(~indices) = posCORR.M(:,mf)*1000;
    
    % explicitly assign colors to vertices, using HCP colormap 
    clrs = connMapVibModes2CmapRainbow_v2(min(plot_data_mf), max(plot_data_mf), 1000);
    clr_IDs = round( (plot_data_mf-min(plot_data_mf)) /max(plot_data_mf-min(plot_data_mf)) *999)+1;
    plot_colors = clrs(clr_IDs,:);
    
    plot_data_mf_w_borders = plot_HCP_boundaries(plot_colors,'all','white');
    
    if ismember(mf,[2,3,4,5,8,10,12]) % plot borders of manually selected regions will reproduce what is shown in Fig.2 
        
        %plot_data_mf_w_borders = plot_HCP_boundaries(plot_colors,'','white');
        plot_data_mf_w_borders = plot_colors;
        if  ~strcmp(border_color,'none')
            % add borders for each dimension (some fct harmonics called more than once)
            if ismember(mf,[2,3]) % visual and somatosensory/motor networks
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'VIS',border_color);
            end
            if mf==3
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'SSM_yeo_som',border_color);
            end
            if mf==2
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'SSM_yeo',border_color);
            end
            if ismember(mf,[3,4,5,8,12]) % somatotopic regions
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'somatotopy',border_color);
            end
            if ismember(mf,[8,12])
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'24dd',border_color);
            end
            if ismember(mf,[5,8]) % early visual
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'earlyVIS',border_color);
            end
            if ismember(mf,7) % early visual
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'ventral_stream',border_color);
            end
            if mf==10
                [plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_data_mf_w_borders,'DMN_yeo',border_color);
            end
        end
        % plot different views 
        % if flat, only left and right
        if strcmp(surf_type,'flat')
            subtightplot(1,nsub,1,[gapleft gapright], margin_height_flat, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            view([5 90])
            material('dull')
            
            subtightplot(1,nsub,2,[gapleft gapright], margin_height_flat, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            view([-10 90])
            material('dull')
        
        % if surface, 6 panels (2xmedial, 2xlateral, top, bottom)
        elseif strcmp(surf_type,'inflated_MSMAll')
            % left medial
            subtightplot(1,nsub,1,[gapleft gapright], margin_height, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            view([90 0])
            material('dull')
            camlight(-190,240)
            
            % right lateral
            subtightplot(1,nsub,2,[gapleft gapright], [0.2 0.2], [0.04 0.02]);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            view([90 0])
            material('dull')
            camlight(145,215)
            
            % both from the top
            subtightplot(1,nsub,3,[gapleft gapright], [0.1 0.1], [0.02 0]);
            options.face_vertex_color = plot_data_mf_w_borders;
            plot_mesh(vertices.all, faces.all, options);
            view([0 90])
            material('dull')
            camlight(0,270)
            
            % both from below
            subtightplot(1,nsub,4,[gapleft gapright], [0.1 0.1], [0 0.04]);
            options.face_vertex_color = plot_data_mf_w_borders;
            plot_mesh(vertices.all, faces.all, options);
            view([180 -90])
            material('dull')
            camlight(180,180)
            
            % left lateral
            subtightplot(1,nsub,5,[gapleft gapright], margin_height, [0.04 0.07]);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            view([270 0])
            material('dull')
            camlight(190,-180)
            
            % right medial
            subtightplot(1,nsub,6,[gapleft gapright], margin_height, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            view([270 0])
            material('dull')
            camlight(-190,-180)
        end
        % plot clusters from kmeans
    elseif ismember(mf,[6,7,9,11]) 
        if mf~=11
            k=2;
        else
            k=3;
        end
        
        c1 = choose_clusters(mf,1);
        if  ~strcmp(border_color,'none')
            % plot one cluster on left medial and right lateral surface
            [plot_data_mf_w_borders,~] = plot_HCP_boundaries_for_clust(plot_colors,mf-1,k,c1,false);
        end
        
        % if flat, only one surface with these cluster borders
        if strcmp(surf_type,'flat')
            subtightplot(1,nsub,1,[gapleft gapright], margin_height_flat, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            view([5 90])
            material('dull')
            % if surface, 1xmedial, 1xlateral with borders of one cluster
        elseif strcmp(surf_type,'inflated_MSMAll')
            % left medial
            subtightplot(1,nsub,1,[gapleft gapright], margin_height, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            camlight;
            view([90 0])
            material('dull')
            camlight(-190,240)
            
            % right lateral
            subtightplot(1,nsub,2,[gapleft gapright], [0.2 0.2], [0.04 0.02]);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            camlight;
            view([90 0])
            material('dull')
            camlight(-190,180)
        end
        
        % plot the other cluster on left lateral and right medial surface
        c2 = choose_clusters(mf,2);
        if ~strcmp(border_color,'none')
            [plot_data_mf_w_borders,~] = plot_HCP_boundaries_for_clust(plot_colors,mf-1,k,c2,false);
        end
        
        % if flat, only one surface with these cluster borders
        if strcmp(surf_type,'flat')
            subtightplot(1,nsub,2,[gapleft gapright], margin_height_flat, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            view([-10 90])
            material('dull')
            % if surface, 1xmedial, 1xlateral with borders of other cluster
        elseif strcmp(surf_type,'inflated_MSMAll')
            % left lateral
            subtightplot(1,nsub,5,[gapleft gapright], margin_height, [0.04 0.07]);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            camlight;
            view([270 0])
            material('dull')
            camlight(-190,180)
            
            % right medial
            subtightplot(1,nsub,6,[gapleft gapright], margin_height, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            camlight;
            view([270 0])
            material('dull')
            camlight(-190,180)
        end
        
        if strcmp(surf_type,'inflated_MSMAll')
            % in the middle plots with both hemispheres, display one cluster on one side, one on the other
            [plot_data_mf_w_borders,~] = plot_HCP_boundaries_for_clust(plot_colors,mf-1,k,[c1,c2],true);
            % both from the top
            subtightplot(1,nsub,3,[gapleft gapright], [0.1 0.1], [0.02 0]);
            options.face_vertex_color = plot_data_mf_w_borders;
            plot_mesh(vertices.all, faces.all, options);
            camlight;
            view([0 90])
            material('dull')
            camlight(0,0)
            
            % both from below
            [plot_data_mf_w_borders,~] = plot_HCP_boundaries_for_clust(plot_colors,mf-1,k,[c2,c1],true);
            subtightplot(1,nsub,4,[gapleft gapright], [0.1 0.1], [0 0.04]);
            options.face_vertex_color = plot_data_mf_w_borders;
            plot_mesh(vertices.all, faces.all, options);
            camlight;
            view([180 -90])
            material('dull')
            camlight(0,0)
        end
        
    end
    if plot_colorbar
        subtightplot(1,nsub,7,[0.04 0.04], [0.025 0], [0.05 0.01]);
        colormap(clrs)
        axis off
        h=colorbar();
        set(h,'XTick',[0,abs(min(plot_data_mf))/(abs(min(plot_data_mf))+abs(max(plot_data_mf))),1])
        % set labels
        xlabs = [min(plot_data_mf),0,max(plot_data_mf)];
        xlabs = round(xlabs);
        set(h,'XTickLabel',xlabs)
        set(h, 'Position', [0.8542    0.2373    0.0179    0.5763]);
        if strcmp(surf_type,'flat')
            set(h,'AxisLocation','out')
        else
            set(h,'AxisLocation','in')
        end
    end
  
    set(gcf,'Color','white')
    f.Units='centimeters';    
    if strcmp(surf_type,'inflated_MSMAll')
        f.Position = [26.4319   19.0235   31.8823    6.8263];
    else
        f.Position = [26.4319   20.2406   20.3465    6.8263];
    end
end





