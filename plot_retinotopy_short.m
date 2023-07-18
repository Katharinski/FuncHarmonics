% this script correlates the functional harmonics with retinotopic
% gradients in visual areas V1-V4, and produces the plots shown in figure
% 3b,c
% code by Katharina Glomb, 2019-2021
% katharina.glomb@gmail.com 

clearvars

homedir = getenv('HOME');
repodatadir = fullfile(homedir,'repos','FuncHarmonics','data');

load('Ind_S900','indices')
nvox = sum(~indices); % #non-medial wall voxels

% extract relevant data from "allresults"
% from documentation: "'allresults' is 91282 grayordinates x 6 quantities x
% 184 datasets x 3 model fits with the full set of pRF analysis results."
% 181 individual subjects + 3 group-average subjects = 184 datasets. The 
% [group] subjects are placed in the 182nd, 183rd, and 184th slots,
% respectively." 
% --> "subject" 182 used all subjects, so we select only this
% "We provide results obtained for each fit (in [this] order [...]: 1=all, 
% 2=first-half-of-each-run, 3=second-half-of-each-run)" 
% --> we select only "all" 
if ~exist('prfresults_group.mat','file')
    load('prfresults.mat','allresults','quants')
    %groupresults = squeeze(allresults(~indices,:,184,1));
    groupresults = squeeze(allresults(1:sum(~indices),:,184,1));
    comment = 'Obtained from prfresults_group.mat, https://osf.io/bw9ec/, groupresults=squeeze(allresults(1:sum(~indices),:,182,1)).';
    save(fullfile(repodatadir,'HCP_derived','prfresults_group.mat'),'groupresults','quants','comment')
else
    load('prfresults_group.mat','groupresults','quants')
end

path_to_surfaces = fullfile(repodatadir,'HCP');
subject_name = 'S900';
% surf_type = 'inflated_MSMAll';
surf_type = 'flat';
try
    [vertices, faces] = connRSMreadGII(path_to_surfaces, subject_name, surf_type); % get surfaces
catch
    error('Error reading surface. Perhaps gifti() function is not working. Try using gifti toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
end

% load fct. harmonics  
load('HCP_S900_CORR_manifold_knn300.mat','posCORR') % posCORR.M contains fct harmonics

%% 
% RETINOTOPY 
% eccentricity: 
ret_data_ecc = groupresults(:,strcmp(quants,'ecc'));
% angle: 
ret_data_ang = groupresults(:,strcmp(quants,'ang'));

%% plot flat maps for display next to retinotopies 
save_now = false;
plot_cbar = true;
mf = 9;
% get retinotopy data and a manifold data
fs = 15; % fontsize
twocolwidth = 18.3; % width of two col figures in cm

plot_data_mf = zeros(length(indices),1);
plot_data_mf(~indices) = posCORR.M(:,mf);
clrs_mf = connMapVibModes2CmapRainbow_v2(min(plot_data_mf), max(plot_data_mf), 1000);
clr_IDs = round( (plot_data_mf-min(plot_data_mf)) /max(plot_data_mf-min(plot_data_mf)) *999)+1;
plot_colors = clrs_mf(clr_IDs,:);
plot_colors(indices,:) = repmat([0.7 0.7 0.7],[sum(indices),1]); % CC is grey
[plot_data_mf_w_borders,~,~] = plot_HCP_boundaries(plot_colors,'earlyVIS','white');

f = figure;
f.Units = 'centimeters';
if plot_cbar
    subplot(211)
end
plot_data_mf_left = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
options.face_vertex_color = plot_data_mf_left;
plot_mesh(vertices.left, faces.left, options);

%view(15,90)
%camlight(110,-90)
material('dull')

if plot_cbar
    subplot(212)
    axis off
    h=colorbar('Location','North','AxisLocation','In');
    h.Position = [0.5000    0.5108    0.3312    0.0188];
    colormap(h,clrs_mf)
    %ylabel(h,{'Functional';'harmonic'})
    ylabel(h,sprintf('$\\psi_{%i}$',mf-1),'interpreter','latex')
    % set labels positions
    %set(h,'XTick',[0,abs(min(plot_data_mf))/(abs(min(plot_data_mf))+abs(max(plot_data_mf))),1])
    set(h,'XTick',[0,1])
    % set labels
    %xlabs = [min(plot_data_mf),0,max(plot_data_mf)];
    xlabs = [min(plot_data_mf),max(plot_data_mf)];
    xlabs = round(xlabs*1000)/1000;
    set(h,'XTickLabel',xlabs)
    set(h,'FontSize',fs)
    f.Position = [11.4375    5.8958    twocolwidth/4    9.2];
else
    f.Position = [11.4375    5.8958  twocolwidth/4    4];
end

set(gcf,'Color','white')

%% polar plot 
% get IDs for visual areas from atlas
load('atlas','glasser2016')
atlas_cortex = glasser2016(1:nvox);
early_vis = {'V1','V2','V3','V4'};
vis_IDs = [1,4,5,6];
atlas_vis_IDs = false(nvox,length(vis_IDs),2);
for v=1:length(vis_IDs)
    atlas_vis_IDs_all = atlas_cortex==vis_IDs(v);
    atlas_vis_IDs(1:nvox/2,v,1) = atlas_vis_IDs_all(1:nvox/2);
    atlas_vis_IDs(nvox/2+1:end,v,2) = atlas_vis_IDs_all(nvox/2+1:end);
end

% map manifold values to colors
% record correlations
cmp = connMapVibModes2CmapRainbow_v2(min(posCORR.M(:,mf)), max(posCORR.M(:,mf)), 990);
clrs = round(((posCORR.M(:,mf)-min(posCORR.M(:,mf)))/max(posCORR.M(:,mf)-min(posCORR.M(:,mf))))*989)+1;
f = figure;
f.Units = 'centimeters';
f.PaperPositionMode = 'auto';
% in order for transparency to work, all points have to plotted in the same
% round
angle_pi = zeros(sum(atlas_vis_IDs(:)),1);
x = zeros(size(angle_pi));
IDs = zeros(size(angle_pi));
counter = 1;
for hem=[1,2]
    for id=1:length(vis_IDs)
        angle_pi_loc = ret_data_ang(atlas_vis_IDs(:,id,hem))/180*pi; % angle in rad
        angle_pi(counter:counter+length(angle_pi_loc)-1) = angle_pi_loc; 
        x_loc = ret_data_ecc(atlas_vis_IDs(:,id,hem)); % eccentricity
        x(counter:counter+length(angle_pi_loc)-1) = x_loc;
        IDs_loc = find(atlas_vis_IDs(:,id,hem));
        IDs(counter:counter+length(angle_pi_loc)-1) = IDs_loc;
        counter = counter+length(angle_pi_loc);
    end
end

for m=1:length(IDs)
    pl = polarplot(angle_pi(m)+rand(1)*0.001,x(m)+rand(1)*0.01,'o', ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cmp(clrs(IDs(m)),:),...
        'MarkerSize',4);
    hold on
end

rlim([0 8])
rticks([0 2 4 6 8])
rticklabels({'','','','',''})
thetaticks(0:90:315)
thetaticklabels({['0',char(176)],['90',char(176)],['180',char(176)],['270',char(176)]})
set(gca,'FontName','FreeSans')
set(gcf,'Color','white')
set(gca,'FontSize',26)
f.Position = [11.4375    5.8958    twocolwidth/6    twocolwidth/6];
f.Position = [25 25 8.9/2 8.9/2];

%% correlations between retinotopic data and manifolds
corrs_ecc = zeros(11,1);
pvals_ecc = zeros(11,1);
corrs_ang = zeros(11,1);
pvals_ang = zeros(11,1);
for mf=2:101%12
    [cc,p] = corrcoef(posCORR.M(IDs,mf),angle_pi);
    corrs_ang(mf-1) = cc(2);
    pvals_ang(mf-1) = p(2);
    [cc,p] = corrcoef(posCORR.M(IDs,mf),x);
    corrs_ecc(mf-1) = cc(2);
    pvals_ecc(mf-1) = p(2);
end


