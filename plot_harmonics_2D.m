%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for plotting two functional harmonics against each other, as 
% shown in the manuscript, Figure 3a, and SI Figure ...
% 
% Copyright (c) 2004 Selen Atasoy
% Author: Selen Atasoy (selenatasoy@gmail.com)
% Date: 21/05/2018
% 
% modified by Katharina Glomb, 2019
% katharina.glomb@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars

homedir = getenv('HOME');
repodatadir = fullfile(homedir,'repos','FuncHarmonics','data');

% load the manifold
load('HCP_S900_CORR_manifold_knn300.mat') % posCORR.M contains the functional harmonics

% load the surfaces
file_path = fullfile(repodatadir,'HCP');
subject_name = 'S900';
surf_type = 'inflated_MSMAll';
try
    [vertices, faces] = connRSMreadGII(file_path, subject_name, surf_type); % get surfaces
catch 
    error('Error reading surface. Perhaps gifti() function is not working. Try using gifti toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
end
% load the indices

load('Ind_S900.mat') % medial wall indices 

% load the labels
if ~exist('HCP_plot_labels.mat','file')
    error('File with surface labels of HCP parcellation does not exist. Get it by running get_surface_labels.m')
else
    load('HCP_plot_labels','surface_labels')
end
HCP.all = surface_labels;
% create the manifold
M = zeros(length(vertices.all), size(posCORR.M,2));
M(~(indices), :) = posCORR.M;

% create the colormap
load('atlas.mat')
c_map= [glasser2016cmap; glasser2016cmap(2:end,:)];

HCP_names = cell(1,361);
HCP_names(1:181) = glasser2016labels;
HCP_names(182:361) = glasser2016labels(2:end);

%% create the plot
% plot two particular ones
twocolwidth = 18.3; % width of two col figures in cm
polar = 0;
% in the paper, we used 
% Fig 3: 3 and 11 
d1 = 3;
d2 = 11;
h = connRSMplotManifold(M, [d1+1, d2+1], HCP.all, c_map, polar);
axis square
axis tight
ax = h.Children(1);
%xlabel({'Functional'; sprintf('harmonic %i',d1-1)})
xlabel(sprintf('$\\psi_{%i}$',d1),'interpreter','latex')
ylabel(sprintf('$\\psi_{%i}$',d2),'interpreter','latex')
%zlabel(sprintf('Harmonic %i',d3-1))
%h.Units = 'centimeters';
%h.Position = [25 25 8.9/2 8.9/2];
h.Units = 'centimeters'; % 2 column figwidth: 183mm, /4=4.5cm
h.Position = [11.4375    5.8958    twocolwidth/5    twocolwidth/5];
set(gca,'LineWidth',2)

opt.FontSize = 8;
opt.LegendBox = 'off';
opt.YMinorTick = 'off';
opt.XMinorTick = 'off';
h = fancy_figure(h, opt);
box off

set(gca,'FontName','FreeSans')





