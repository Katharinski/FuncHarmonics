% compute modified silhouette values for testing correspondence between
% functional harmonics (+ control basis function sets) and Glasser
% parcellation
%
% Katharina Glomb 
% last modified June 23, 2021
% katharina.glomb@gmail.com 

clearvars

plot_demo = false; % plot a demo of the spherical rotations

% parameters/switches
% which set of basis functions to test
test_which = 'harmonics'; % harmonics, eigvecs, FC_eigvecs, PCs, ICs15
suffix = ''; % optional suffix for testing purposes 
path_to_surfaces = [pwd,'/data'];

if ~exist('HCP_plot_labels.mat','file')
    error('File with surface labels of HCP parcellation does not exist. Get it by running get_surface_labels.m')
else
    load('HCP_plot_labels','surface_labels')
end

%% set up file names 
switch test_which 
    case 'harmonics'
        load('HCP_S900_CORR_manifold_knn300.mat') 
        use_dims = (2:12); % discard constant harmonic
    case 'eigvecs'
        load('HCP_S900_CORR_eigvecs_knn300.mat')
        use_dims = (1:11);
    case 'FC_eigvecs'
        load('eigs_denseFC.mat','V')
        posCORR.M = V;
        use_dims = (1:11);
    case 'PCs'
        load('HCP_group_eigenmaps','posCORR')
        use_dims = (1:11);
    case 'ICs15'
        load('ICs15.mat','ICs')
        posCORR.M = ICs;
        use_dims = (1:11);
end

save_fname_emp = ['silh_coeffs_',test_which,suffix,'.mat'];
save_fname_rand = ['silh_coeffs_spherical_',test_which,suffix,'.mat'];
        
%% compute non-randomized silh vals

if ~exist(save_fname_emp,'file')
    fprintf('Computing empirical silhouette coefficients...\n')
    % for each dimension separately
    [SC_all, M_within, M_between] = compute_silh_vals(posCORR.M(:,use_dims), surface_labels);
    
    % all dimensions together
    [SC_ddim, ~,~] = compute_silh_vals(posCORR.M(:,use_dims),surface_labels,true);
    save(save_fname_emp,'SC_all','M_within','M_between','SC_ddim')
    fprintf('Done!\n')
else
    load(save_fname_emp,'SC_all','M_within','M_between')
end

%% compute silh vals of rotations
alpha = 0.05;
if ~exist(save_fname_rand,'file')
    fprintf('Computing rotated silhouette coefficients...\n')
    % determine number of randomizations
    nrand = (1/alpha) * length(use_dims);
    % inititalize arrays
    parcels = unique(surface_labels);
    if parcels(1)==0 % zero label will be considered an "unknown" region, to be ignored
        parcels = parcels(2:end);
    end
    nparc = length(parcels);
    % each dim separately
    SC_all_rand = zeros(nparc,length(use_dims),nrand); 
    M_within_rand = zeros(nparc,length(use_dims),nrand); 
    M_between_rand = zeros(nparc,length(use_dims),nrand); 
    % all dimensions together
    SC_ddim_rand = zeros(nparc,nrand); 
    % load or create rotations
    if ~exist('rotations_harmonics.mat','file')
        rotations = create_rotations(nrand,path_to_surfaces); 
        save('rotations_harmonics.mat','rotations')
    else
        load('rotations_harmonics','rotations')
        if size(rotations,2)<nrand
            warning('Number of pre-computed rotations is smaller than indicated number of randomizations.')
        end
    end
    
    nV = size(rotations,1); % number of vertices including CC
    load('Ind_S900','indices')
    for n=1:nrand
        fprintf('Randomization %i of %i...\n',n,nrand)
        dim_counter = 0;        
        M_rand = zeros(size(posCORR.M,1),length(use_dims)); % save rotated data, but w/o CC vertices
        % indices for rotations
        Il = rotations(1:nV/2,n);
        Ir = rotations((nV/2+1):end,n);
        % apply to data
        for mf=use_dims
            dim_counter = dim_counter+1;
            % CC needed for this to work (otherwise spherical proj not valid)
            dataLR = nan(nV,1); 
            dataLR(~indices) = posCORR.M(:,mf); % CC will remain NaN
            dataL = dataLR(1:nV/2);
            dataR = dataLR((nV/2+1):end);
            dataLrot = dataL(Il);
            dataRrot = dataR(Ir); 
            
            % compute silhouette values
            % remove CC for this computation
            dataLRrot = [dataLrot;dataRrot];
            dataLRrot = dataLRrot(~indices,:); % remove CC for computing silhouette values 
            M_rand(:,dim_counter) = dataLRrot; % save for using multidim distance after this loop ends
        end
        % compute the silhouette values using rotated maps
        % for each dimension separately
        [SC_all_rand(:,:,n), M_within_rand(:,:,n), M_between_rand(:,:,n)] = compute_silh_vals(M_rand, surface_labels);
        % in d-dimensional space
        [SC_ddim_rand(:,n), ~,~] = compute_silh_vals(M_rand,surface_labels,true);
    end
    save(save_fname_rand,'SC_all_rand','M_within_rand','M_between_rand','SC_ddim_rand')
    fprintf('Done!\n')
else
    load(save_fname_rand)
    nrand = size(SC_all_rand,3);
end

%% demo 
% plot before and after rotating 
% load regular surfaces
if plot_demo==true
    if strcmp(getenv('USER'),'katha')
        working_dir = '/home/katha/';
    elseif strcmp(getenv('USER'),'localadmin')
        working_dir = '/home/localadmin/';
    elseif strcmp(getenv('USER'),'katharina')
        working_dir = '/mnt/data/repos/';
    end
    
    subject_name = 'S900';
    surf_type = 'inflated_MSMAll';
    [vertices, faces] = connRSMreadGII(path_to_surfaces, subject_name, surf_type); % get surfaces
    
    mf = 2; % choose 2:12
    % last rotation will be used
    dataLR = nan(nV,1);
    dataLR(~indices) = posCORR.M(:,mf); % CC will remain NaN
    dataL = dataLR(1:nV/2);
    dataR = dataLR((nV/2+1):end);
    dataLrot = dataL(Il);
    dataRrot = dataR(Ir);
    
    f1 = plot_HCP_on_surface([dataL;dataR],vertices,faces);
    pos = f1.Position;
    f1.Position = [pos(1) pos(2) 20 4];
    f2 = plot_HCP_on_surface([dataLrot;dataRrot],vertices,faces);
    pos = f2.Position;
    f2.Position = [pos(1) pos(2) 20 4];
end

%% plot real and rotated silhouette values 
av_SC_emp = squeeze(mean(SC_all));
av_SC_rand = squeeze(nanmean(SC_all_rand));

% clrs = [230,25,75;128,128,0;0,128,128;170,110,40]/255; % red, olive, teal,brown
clrs = [230,25,75;181,82,106;118,81,89;0,128,128;0,0,0]/255; % pink, another pink, another pink, teal 
grayclr = [0.65, 0.65, 0.65]; 

gapvert = 0.01;
gaphorz = 0.07;
margin_height = [0.01 0.05]; % lower, upper
margin_width = [0.2 0.05]; % left, right

switch test_which 
    case 'harmonics'
        xlab = {'Functional';'harmonics'};
        leglab = 'Functional harmonics';
        clr = clrs(1,:); 
    case 'eigvecs'
        xlab = {'Adjacency'; 'eigenvectors'};
        leglab = 'Adj. eigenvectors';
        clr = clrs(2,:); 
    case 'FC_eigvecs'
        xlab = {'FC eigenvectors';'  '};
        leglab = 'FC eigenvectors';
        clr = clrs(3,:); 
    case 'PCs'
       xlab = {'Principal';'components'};
       leglab = 'Principal components';
       clr = clrs(4,:); 
    case 'ICs15'
        xlab = {'Independent'; 'components'};
        leglab = 'Ind. components';
        clr = 'k';
end

plot_use_dims = 1:11; % though we skip the first harmonic, it's constant, so we call it the "0th harmonic" and count only non-constant ones when plotting 

twocolwidth = 18.3; % cm
f = figure;
subtightplot(7,1,[1:5],[gapvert gaphorz], margin_height, margin_width);
h2=plot(plot_use_dims,av_SC_rand,'x','Color',grayclr);
hold on
h1=plot(plot_use_dims,av_SC_emp,'o','Color',clr,'LineWidth',2);
set(gca,'FontName','FreeSans')
set(gca,'XLim',[0.9 length(plot_use_dims)+0.1])
set(gca,'YLim',[0.4 0.9])
outpos = get(gca,'OuterPosition'); 
set(gca,'OuterPosition',[outpos(1) 0.3 outpos(3) 0.7])
xlabel(xlab)
ylabel('Mod. silhouette value')
set(gca,'FontSize',8)
f.Units='centimeters';
pos = f.Position;
set(gcf,'Resize','off')
f.Position = [pos(1) pos(2) twocolwidth/5 5.5];
%f.Position = [pos(1) pos(2) 7.7 6];

box off
grid on
% legend in extra subplot
subtightplot(7,1,7,[gapvert gaphorz], margin_height, margin_width);
plot(1,1,'o','Color',clr,'LineWidth',2);
hold on
plot(1,1,'x','Color',grayclr);
set(gca,'XLim',[0,0.5]); % plotted point will be invisible
axis off
[l,icons,~] = legend({leglab,'Rotated'},'Location','Best','EdgeColor','white','FontName','FreeSans');
l.Position = [-0.045    0.1026    1.2000    0.1345];
icons(1).FontName = 'FreeSans';
icons(2).FontName = 'FreeSans';
% adjust distance between icon and text
icons(1).Position = [0.2 icons(1).Position(2) 0];
icons(2).Position = [0.2 icons(2).Position(2) 0];
icons(1).FontSize = 8;
icons(2).FontSize = 8;
set(gcf,'Color','white')
f.Position = [pos(1) pos(2) twocolwidth/5 5.5];

save_now = false;

if save_now
    export_fig(['average_SCs_',test_which,'_small',suffix,'_goodfont.png'],'-r300','-nocrop','-opengl');
    export_fig(['average_SCs_',test_which,'_small',suffix,'_goodcircles.png'],'-r300','-nocrop','-painters');    
end

%% plot everything together 
test_all = cell(5,1); 
test_all{1} = 'harmonics';
test_all{2} = 'eigvecs';
test_all{3} = 'FC_eigvecs';
test_all{4} = 'PCs';
test_all{5} = 'ICs15';



%% Wilcoxon ranksum test 
% can't use signed rank bc components can't be matched 
% Benjamini-Hochberg to correct for multiple comparisons
B = 5;
ncomp = nchoosek(B,2); % pairs of bases
    
all_p = zeros(ncomp,1); 
% record which bases are compared and which one won (z-value)
all_zval = zeros(ncomp,3);

counter = 1; 
for b1=1:B
    
    load_fname1 = ['silh_coeffs_',test_all{b1},'.mat'];
    load(load_fname1,'SC_all')
    set1 = squeeze(mean(SC_all))';
    
    for b2=(b1+1):B
        all_zval(counter,1) = b1; 
        all_zval(counter,2) = b2;
        load_fname2 = ['silh_coeffs_',test_all{b2},'.mat'];
        load(load_fname2,'SC_all')
        set2 = squeeze(mean(SC_all))';
    
        [all_p(counter),h,stats] = ranksum(set1,set2);
        % compute signrank myself since Matlab returns something funny
        all_zval(counter,3) = sign(stats.zval);
        counter = counter+1; 
    end
end

crit_p_BH = (([1:ncomp]./ncomp)*0.05)';
[sortP,orderP] = sort(all_p);
sort_all_zval = all_zval(orderP,:);
sig = sortP<crit_p_BH;

all_sig_diffs = zeros(sum(sig),3);
all_sig_diffs(:,[1,2]) = sort_all_zval(sig,[1,2]);
all_sig_diffs(:,3) = sortP(sig);

short_model_names = test_all;

if sum(sig)>0
    fprintf('significant results: \n')
    sig = find(sig);
    for sg=1:length(sig)
        if sort_all_zval(sg,3) <0
            fprintf('%s is better than %s\n',short_model_names{sort_all_zval(sg,2)},short_model_names{sort_all_zval(sg,1)})
        elseif sort_all_zval(sg,3) >0
            fprintf('%s is better than %s\n',short_model_names{sort_all_zval(sg,1)},short_model_names{sort_all_zval(sg,2)})
        end
    end
    fprintf('-------\n')
else
    fprintf('No significant results :(\n')
end

%% plot
% clrs = [230,25,75;128,128,0;0,128,128;170,110,40;0,0,0]/255;

% twocolwidth = 18.3; % cm
f2 = figure;
subtightplot(7,1,[1:5],[gapvert gaphorz], margin_height, margin_width);
hold on
for bs=1:5
    load_fname = ['silh_coeffs_',test_all{bs},'.mat'];
    load(load_fname,'SC_all')
    av_SC_emp = squeeze(mean(SC_all));
    plot(bs,av_SC_emp,'o','Color',clrs(bs,:),'LineWidth',2)
end
set(gca,'FontName','FreeSans')
set(gca,'XLim',[0.5 B+0.5])
set(gca,'YLim',[0.4 0.9])
set(gca,'XTick',1:5)
outpos = get(gca,'OuterPosition'); 
set(gca,'OuterPosition',[outpos(1) 0.3 outpos(3) 0.7])
xlabel(xlab)
ylabel('Mod. silhouette value')
set(gca,'FontSize',8)
f2.Units='centimeters';
set(gcf,'Resize','off')
f2.Position = [pos(1) pos(2) twocolwidth/5 5.5];

box off
grid on
% legend in extra subplot
subtightplot(7,1,7,[gapvert gaphorz], margin_height, margin_width);
plot(1,1,'o','Color',clr,'LineWidth',2);
hold on
plot(1,1,'x','Color',grayclr);
set(gca,'XLim',[0,0.5]); % plotted point will be invisible
axis off
[l,icons,~] = legend({leglab,'Rotated'},'Location','Best','EdgeColor','white','FontName','FreeSans');
l.Position = [-0.045    0.1026    1.2000    0.1345];
icons(1).FontName = 'FreeSans';
icons(2).FontName = 'FreeSans';
% adjust distance between icon and text
icons(1).Position = [0.2 icons(1).Position(2) 0];
icons(2).Position = [0.2 icons(2).Position(2) 0];
icons(1).FontSize = 8;
icons(2).FontSize = 8;
set(gcf,'Color','white')
% f.Position = [pos(1) pos(2) twocolwidth/5 5.5];

% -nocrop necessary to keep size
% super cool bug (loads of fun): when using opengl, the circles look awful,
% but when using painters, the font in the legend isn't correct, so I need
% to combine both afterwards! SSOOOO FFFUUUUUNNNNNN

save_now = false;

if save_now
    export_fig('average_SCs_all_small_goodfont.png','-r300','-nocrop','-opengl');
    export_fig('average_SCs_all_small_goodcircles.png','-r300','-nocrop','-painters');    
end

%%

f3 = figure; 
lbls = test_all; 
lbls{3} = 'FC eigves';
lbls{5} = 'ICs'; 

hold on
for bs=1:5
    load_fname = ['silh_coeffs_',test_all{bs},'.mat'];
    load(load_fname,'SC_all')
    av_SC_emp = squeeze(mean(SC_all));
    plot(bs,av_SC_emp,'o','Color',clrs(bs,:),'LineWidth',2)
end

f3.Units='centimeters';

set(gca,'FontName','FreeSans')
set(gca,'XLim',[0.5 B+0.5])
ylim = get(gca,'YLim'); 
set(gca,'YLim',[0.4 ylim(2)+0.19])
outpos = get(gca,'OuterPosition'); 
set(gca,'OuterPosition',[outpos(1) 0.3 outpos(3) 0.7])
ylabel('Mod. silhouette value')
set(gcf,'Resize','off')
f3.Position = [pos(1) pos(2) twocolwidth/5 6.5];
set(gca,'XTick',1:5)
set(gca,'XTickLabel',lbls)
set(gca,'FontSize',8)
RotateXLabel(90,lbls);
% grid on

% add significance stars 
% sort so it looks good 
[~,sortID] = sort(all_sig_diffs(:,1)); 
all_sig_diffs_sort = all_sig_diffs(sortID,:); 
all_sig_diffs_sort2 = zeros(size(all_sig_diffs)); 
sig1 = unique(all_sig_diffs(:,1));
counter = 1; 
for en=1:length(sig1)
    nen = sum(all_sig_diffs(:,1)==sig1(en)); 
    this_en = all_sig_diffs(all_sig_diffs(:,1)==sig1(en),:); 
    [~,sortID2] = sort(this_en(:,2)); 
    all_sig_diffs_sort2(counter:(counter+nen-1),:) = this_en(sortID2,:); 
    counter = counter+nen; 
end
all_sig_diffs = all_sig_diffs_sort2;
ini_shift = 0.2; 
len_bars = 0.02; 
shift = 0.03; 
for sd=1:length(all_sig_diffs)
    plot([all_sig_diffs(sd,1),all_sig_diffs(sd,2)],[ylim(2)+ini_shift-(sd*shift),ylim(2)+ini_shift-(sd*shift)],'k')
    plot([all_sig_diffs(sd,1),all_sig_diffs(sd,1)],[ylim(2)+ini_shift-(sd*shift)-len_bars,ylim(2)+ini_shift-(sd*shift)],'k')
    plot([all_sig_diffs(sd,2),all_sig_diffs(sd,2)],[ylim(2)+ini_shift-(sd*shift)-len_bars,ylim(2)+ini_shift-(sd*shift)],'k')
end

% pos = f2.Position;
set(gcf,'Color','white')
% f3.Position = [pos(1) pos(2) twocolwidth/5 5.5];

save_now = false;

if save_now
    export_fig('average_SCs_all_small_xlab_goodfont.png','-r300','-nocrop','-opengl');
    export_fig('average_SCs_all_small_xlab_goodcircles.png','-r300','-nocrop','-painters');    
end

