% script for plotting auditory hierarchy as shown in Figure S... 
% code by Katharina Glomb, 2019 
% katharina.glomb@gmail.com 

clearvars
% grouping of cortical areas according to Glasser et al. 2016 supplementary
% neuroanatomical results 
early_aud = {'A1','MBelt','LBelt','PBelt','RI'};
early_aud_IDs = [24,173,174,124,104];
ass_aud = {'A4','A5','STSdp','STSda','STSvp','STSva','TA2','STGa'};
ass_aud_IDs = [175,125,129,128,130,176,107,123];

all_AUD_IDs = [early_aud_IDs,ass_aud_IDs];

% load the functional harmonics 
load('HCP_S900_CORR_manifold_knn300.mat','posCORR') % posCORR.M contains the fct. harmonics 
load('Ind_S900.mat')

if ~exist('HCP_plot_labels.mat','file')
    error('File with surface labels of HCP parcellation does not exist. Get it by running get_surface_labels.m')
else
    load('HCP_plot_labels','surface_labels')
end

load('atlas','glasser2016cmap')
mf = 11;
M = zeros(64984,1);
mf_vals = posCORR.M(:,mf);
M(~indices) = mf_vals;
mingrad = min(M(:));
maxgrad = max(M(:));

ffs = 8;
twocolwidth = 18.3; % width of two col figures in cm

%% extract "auditoryness" for all regions based on HCP colormap 
% as described in Glasser et al. 2016, caption of Figure 3: "Colours 
% indicate the extent to which the areas are associated in the resting 
% state with auditory (red) [...] groups of areas."

glasser2016cmap = glasser2016cmap(2:end,:); % remove "unknown" parcel
ncl = size(glasser2016cmap,1);

red = [1 0 0];
reddists = zeros(ncl,1); % auditory = red

for cl=1:ncl
    reddists(cl) =  norm(glasser2016cmap(cl,:)-red);
end
%% plot and save figure 
save_now = false;
f = figure;
hold on
meansparcels = zeros(ncl,1);
for r=1:ncl
    this_parcel = M(surface_labels==r);
    % plot all parcels 
    %plot(repmat(glasser2016cmap(r,1),[length(this_parcel),1]),this_parcel,'.','Color',glasser2016cmap(r,:))
    % plot only the mean of each parcel 
    meansparcels(r) = mean(this_parcel);
    plot(glasser2016cmap(r,1),mean(this_parcel),'.','Color',glasser2016cmap(r,:),'MarkerSize',10);%,'LineWidth',2)
end
%cc = corrcoef(glasser2016cmap(:,1)./sum(glasser2016cmap,2),meansparcels);
[cc,p] = corrcoef(glasser2016cmap(:,1),meansparcels);
xlabel({'Association with';'auditory areas'})
%xlabel(sprintf('$\\psi_%i$',d1-1),'interpreter','latex')
%ylabel({'Parcellated';'$\psi_{10}$'},'interpreter','latex')
ylabel('Parcellated')
text(0.1,0.01,sprintf('r=%.2f,p=%.0d',cc(2),p(2)),'FontSize',ffs,'FontName','FreeSans')
%title(mf-1)

set(gcf,'Color','white')
opt.FontSize = ffs;
opt.LegendBox = 'off';
set(gca,'LineWidth',1)
%h1 = fancy_figure(f1, opt);
%box off
%h2 = fancy_figure(f2,opt);
box off
set(gca,'FontSize',ffs)
set(gca,'FontName','FreeSans')
f.Units = 'centimeters';
f.Position = [11.4375    5.8958    twocolwidth/5    twocolwidth/5]; % h.Position = [11.4375    5.8958    twocolwidth/5    twocolwidth/5];

if save_now
    export_fig('auditoryness_harmonic10.png','-r300');
end