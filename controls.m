% script for testing alternative dimensionality reduction techniques in
% terms of their reconstruction performance 
% The comparison is always between harmonics and a control basis set,
% or the basis set (control basis set or harmonics) and rotations thereof
%
% Katharina Glomb 
% last modified June 25, 2021
% katharina.glomb@gmail.com 

clearvars

%% load control basis set 
test_which = 'harmonics'; % harmonics, eigvecs, PCs, rotations, FC_eigvecs
use_rotations = false; 
clrs = [230,25,75;181,82,106;118,81,89;0,128,128]/255; % pink, another pink, another pink, teal 
grayclr = [0.65, 0.65, 0.65]; 
switch test_which 
    case 'harmonics'        
        load('HCP_S900_CORR_manifold_knn300.mat') % real fct harmonics 
        this_clr = clrs(1,:); 
    case 'eigvecs' % eigenvectors obtained from binarized adjacency matrix
        load('HCP_S900_CORR_eigvecs_knn300.mat','posCORR')
        this_clr = clrs(2,:); 
    case 'PCs' % first 101 PCs provided by HCP
        load('HCP_group_eigenmaps','posCORR')
        this_clr = clrs(4,:); 
    case 'FC_eigvecs'
        load('eigs_denseFC.mat','V')
        posCORR.M = V;
        clearvars V
        this_clr = clrs(3,:); 
    case 'ICs15'
        load('ICs15.mat','ICs')
        posCORR.M = ICs;
        this_clr = 'k';
end

if use_rotations==true
    load('rotations_harmonics','rotations') % indices for rotating fct harmonics - randomization that preserves smoothness
    load('Ind_S900','indices') % CC indices
    V_full = length(indices);
    suffix = '_rotated';
else
    suffix = '';
end

%% load the task maps 
if ~exist('HCP_task_maps.mat','file')
    fname = 'HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s2_MSMAll.dscalar.nii';
    wb_dir = '/home/katharina/snap/workbench/bin_linux64/wb_command';
    try
        ciftiData = ciftiopen(fname,wb_dir);
    catch 
        error('Error reading task maps. Perhaps ciftiopen() function is not working. Try using toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
    end
    taskMaps = ciftiData.cdata;
    save('HCP_task_maps.mat','taskMaps')
else
    load('HCP_task_maps.mat')
end

%% task map reconstruction 
neg_maps_IDs = [12,13,23:26,34,35,50:56,67,68,72,73,79,80,84,85]; % some maps are just inversions of others [?!] - remove them
id_contrasts_IDs = [14,27:30,36,57:62,66,74,78,86];

reduced_flag = false;
if reduced_flag==false
    nT_full = size(taskMaps,2);
    reduced_flag = true;
    taskMaps = taskMaps(:,(~ismember(1:nT_full,[neg_maps_IDs,id_contrasts_IDs])));
    nT = size(taskMaps,2);
end

if strcmp(test_which(1:3),'ICs') % figure out which dimensions to use for the comparison
    use_bf = [1:str2num(test_which(4:end))];
    use_bf_harmonics = [1:21,31:10:101];
    if length(use_bf)<length(use_bf_harmonics)
        matches = ismember(use_bf,use_bf_harmonics);
        if sum(matches)>0
            use_bf = use_bf(matches);
        end
    else
        matches = ismember(use_bf_harmonics,use_bf);
        if sum(matches)>0
            use_bf = use_bf_harmonics(matches);
        end
    end
else
    use_bf = [1:21,31:10:101];
end

%% performance of chosen basis set ('test_which')
if use_rotations==true
    nrot = size(rotations,2);
    nrounds = nrot;
else
    nrounds = 1; 
end

recon_errs_controls_filename = ['RE_corrs','_',test_which,suffix,'.mat'];

if ~exist(recon_errs_controls_filename,'file') 
    % initialize
    if strcmp(test_which(1:3),'ICs')
        spec = zeros(str2double(test_which(4:end)),nT,nrounds);
    else
        spec = zeros(max(use_bf),nT,nrounds);
    end
    RE_all = zeros(length(use_bf),nT,nT,nrounds);
    RE_match = zeros(length(use_bf),nT,nrounds);
    corr_all = zeros(length(use_bf),nT,nT,nrounds);
    corr_match = zeros(length(use_bf),nT,nrounds);
    
    % randomize data if rotations are tested
    for n=1:nrounds
        fprintf('Randomizations, round %i of %i...\n',n,nrounds)
        if use_rotations==true
            dataLR = zeros(V_full,size(posCORR.M,2));
            dataLR(~indices,:) = posCORR.M;
            dataLRrot = dataLR(rotations(:,n),:);
            
            % remove CC for this computation
            this_vectors = dataLRrot(~indices,:);
        else
            this_vectors = posCORR.M; 
        end
        
        [spec(:,:,n),~,RE_all(:,:,:,n),RE_match(:,:,n),corr_all(:,:,:,n),corr_match(:,:,n),~] = compute_recon_errors(this_vectors,taskMaps,use_bf);
    end
    if use_rotations==true
        spec_rand = spec;
        RE_all_rand = RE_all;
        RE_match_rand = RE_match;
        corr_all_rand = corr_all;
        corr_match_rand = corr_match;
        save(recon_errs_controls_filename,'spec_rand','RE_all_rand','RE_match_rand','corr_all_rand','corr_match_rand','use_bf')
    else
        save(recon_errs_controls_filename,'spec','RE_all','RE_match','corr_all','corr_match','use_bf')
    end
else
    load(recon_errs_controls_filename)
end
fprintf('Done!\n')

%% performance of "other" basis set 
% if test_which is harmonics, this cannot be harmonics 

if use_rotations==true % compare to this basis set
    recon_errs_filename = ['RE_corrs','_',test_which,'.mat'];
elseif ~strcmp(test_which,'harmonics')==true % compare to harmonics
    recon_errs_filename = 'RE_corrs_harmonics.mat';
else
    error('Cannot compare harmonics to itself. Choose a different basis set, or use rotations.')
end

if ~exist(recon_errs_filename,'file')
    error('File %i does not exist. Use recon.m to create it.',recon_errs_filename)
else
    use_bf_emp = load(recon_errs_filename,'use_bf'); % first test if correct dimensions were used 
    if ~isequal(use_bf_emp.use_bf,use_bf)
        warning('use_bf not consistent between random and empirical data. Trying to match.')
        if length(use_bf)<length(use_bf_emp.use_bf)
            % which of the use_bf are also in use_bf_emp.use_bf
            matches = ismember(use_bf,use_bf_emp.use_bf);
            % which of the use_bf_emp.use_bf are also in use_bf
            matches2 = ismember(use_bf_emp.use_bf,use_bf);
            if ~strcmp(test_which,'harmonics') && use_rotations==false % in this case, use_bf_emp.use_bf will originate from harmonics and 1st comp is constant - has to be kept bc it's removed later
                matches2 = ismember(use_bf_emp.use_bf-1,use_bf);
                matches2(1) = true; 
            end
            if sum(matches)>0
                use_bf = use_bf(matches);
            end
        else
            matches = ismember(use_bf_emp.use_bf,use_bf);
            if sum(matches)>0
                use_bf = use_bf_emp.use_bf(matches);
            end
        end
    else
        matches = true(length(use_bf),1); 
        matches2 = matches;
    end
    if sum(matches)>0
        if use_rotations==true % variables will have suffix _rand
            load(recon_errs_filename);
            RE_all = RE_all(matches,:,:);
            RE_match = RE_match(matches,:);
            corr_all = corr_all(matches,:,:);
            corr_match = corr_match(matches,:);
        else % variables will have the same name - need to be renamed
            spec_rand = spec;
            RE_all_rand = RE_all(matches,:,:);
            RE_match_rand = RE_match(matches,:);
            corr_all_rand = corr_all(matches,:,:);
            corr_match_rand = corr_match(matches,:);
            load(recon_errs_filename,'RE_all','RE_match','corr_all','corr_match')
            RE_all = RE_all(matches2,:,:);
            RE_match = RE_match(matches2,:);
            corr_all = corr_all(matches2,:,:);
            corr_match = corr_match(matches2,:);
        end
    end    
end

%% plot reconstruction errors for each task type  
save_now = false; 
twocolwidth = 18.3;
fs = 8;
if use_rotations==true && ~exist('RE_match_rand_full','var')
    RE_match_rand_full = RE_match_rand;
    RE_match_rand = min(RE_match_rand,[],3);
end
tasktypes = [ones(30,1);ones(6,1)*2;ones(26,1)*3;ones(6,1)*4;ones(6,1)*5;ones(6,1)*6;ones(6,1)*7];
tasktypes_red = tasktypes(~ismember(1:nT_full,[neg_maps_IDs,id_contrasts_IDs]));
labels = {'Working memory';'Gambling';'Motor';'Language';'Social';'Relational';'Emotion'};
f=figure;
f.Units = 'centimeters';
%clrs = get(gca,'ColorOrder');
for tt=1:7
    subplot(2,4,tt)
    set(gca,'FontName','FreeSans')
    hold on
    h = gobjects(2,1);
    transparent = 1;
    % solid curve: 
    if strcmp(test_which,'ICs15')
        this_use_bf = use_bf_emp.use_bf(matches2); 
        x = this_use_bf(2:end)-1; % non-constant harmonics
        this_task_type = RE_match(2:end,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        try
            H = shadedErrorBar(x',y,errBar,{'Color',clrs(1,:),'LineWidth',2},transparent);
        catch
            error('Error plotting shaded error bar. Perhaps shadedErrorBar function is missing - download at https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar')
        end
    elseif (strcmp(test_which,'harmonics') && (use_rotations==true)) || (~strcmp(test_which,'harmonics') && (use_rotations==false)) % first curve is harmonics in its color
        x = use_bf(2:end)-1; % non-constant harmonics
        this_task_type = RE_match(2:end,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',clrs(1,:),'LineWidth',2},transparent);
    elseif ~strcmp(test_which,'harmonics') && (use_rotations==false) % first curve is other basis set in its color
        x = use_bf(1:end); % first component is not constant and is used 
        this_task_type = RE_match(:,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',clrs(1,:),'LineWidth',2},transparent);
    end
    h(1) = H.mainLine;
    
    % dashed curve 
    if (strcmp(test_which,'harmonics') && (use_rotations==true)) % second curve is harmonics' rotations, in grey 
        x = use_bf(2:end)-1; % non-constant harmonics
        this_task_type = RE_match_rand(2:end,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',grayclr,'LineWidth',2,'LineStyle','--'},transparent);
    elseif ~strcmp(test_which,'harmonics') % second curve is other basis set in its color
        x = use_bf(1:end);
        this_task_type = RE_match_rand(:,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',this_clr,'LineWidth',2,'LineStyle','--'},transparent);
    end
        
    h(2) = H.mainLine;    
%     title(labels{tt})
    if tt==5 || tt==1
        ylabel({'Normalized';'reconstruction error'})
    end
    if strcmp(test_which,'ICs15')
        set(gca,'XLim',[1 15])
    else
        set(gca,'XLim',[1 100])
        set(gca,'XTick',[20:40:100])
        set(gca,'XTickLabel',[20:40:100])
    end
    set(gca,'YLim',[0 1.7])
    xlabel('Components')
    set(gca,'FontSize',fs)
    grid on
end

subplot(2,4,8)
if ~strcmp(test_which,'harmonics')
    plot(1,1,'Color',clrs(1,:),'LineWidth',2);
else
    plot(1,1,'Color',this_clr,'LineWidth',2);
end
hold on
if use_rotations==true && strcmp(test_which,'harmonics')==true
    plot(1,1,'Color',grayclr,'LineWidth',2,'LineStyle','--');
else
    plot(1,1,'Color',this_clr,'LineWidth',2,'LineStyle','--');
end

if use_rotations==false % compare harmonics to controls 
    if strcmp(test_which,'eigvecs')
        [~,icons,~,~] = legend('Fct. harmonics','Adj. eigenvecs','Location','Best');
    elseif strcmp(test_which,'FC_eigvecs')
        [~,icons,~,~] = legend('Funtional harmonics','FC eigenvecs','Location','Best');
    elseif strcmp(test_which,'PCs')
        [~,icons,~,~] = legend('Funtional harmonics','PCs','Location','Best');
    elseif strcmp(test_which(1:3),'ICs')
        [~,icons,~,~] = legend('Funtional harmonics','ICs','Location','Best');
    end
elseif use_rotations==true % compare harmonics and controls to their rotations
    if strcmp(test_which,'harmonics')
        [~,icons,~,~] = legend('Funtional harmonics','Rotations','Location','Best');
    elseif strcmp(test_which,'eigvecs') && use_rotations==true
        [~,icons,~,~] = legend('Adj. eigenvecs','Rotations','Location','Best');
    elseif strcmp(test_which,'FC_eigvecs')
        [~,icons,~,~] = legend('FC eigenvecs','Rotations','Location','Best');
    elseif strcmp(test_which,'PCs')
        [~,icons,~,~] = legend('PCs','Rotations','Location','Best');
    end
end
ylabel({'Normalized';'reconstruction error'})
legend boxoff
set(gca,'FontSize',fs)
% adjust legend line lengths
leglnspos = icons(3).XData;
for ch=(2+1):2:(2*3)
    icons(ch).XData = [leglnspos(2)-0.12 leglnspos(2)];
end
axis off
set(gcf,'Color','white')
f.Position = [26.3790   19.5527   0.75*twocolwidth   8];

if save_now==true && use_rotations==false 
    export_fig(sprintf('SI_fig_RE_tasks_%s.png',test_which),'-r300')
elseif save_now==true && use_rotations==true
    export_fig(sprintf('SI_fig_RE_tasks_%s_rotated.png',test_which),'-r300')
end

%% significance test (permutation test) 
from = 1; 
upto = 15; 

warning('Significance is tested using dimensions %i to %i!',from,upto)

from_ID = find(use_bf==from); 
if isempty(from_ID)
    error('Choose a value for from that occurs in use_bf!')
end

upto_ID = find(use_bf==upto); 
if isempty(upto_ID)
    error('Choose a value for upto that occurs in use_bf!')
end

set1 = RE_match_rand(from_ID:upto_ID,:);
if use_rotations==false || (use_rotations==true && strcmp(test_which,'harmonics')) % RE_match belongs to harmonics: discard constant harmonic
    if upto<=20 % beyond 21, not all dimensions are computed 
        set2 = RE_match(from_ID+1:upto_ID+1,:);
    else
        set2 = RE_match(from_ID:upto_ID,:);
    end
else
    set2 = RE_match(from_ID:upto_ID,:);
end
alpha = 0.05; 
p_corr = (alpha/(7*5)); % 7 task groups, 5 controls 
nperm = 1000; 
pVals = zeros(tt,1); 
realDiff = zeros(tt,1); 

% for each task group
for tt=1:7
    realVals_emp = set2(:,tasktypes_red==tt);
    realVals_rand = set1(:,tasktypes_red==tt); 
    realDiff(tt) = mean(realVals_emp(:)) - mean(realVals_rand(:)); % hypothesis that is tested 
    [pVals(tt), ~] =  checkDiffSignificanceMonteCarlo([realVals_emp(:),realVals_rand(:)], nperm);
end

% output results
nsigresults = 0; 
for tt=1:7
    if (realDiff(tt)<0) && (pVals(tt)<p_corr) % mean empirical RE is significantly smaller than mean control RE
        nsigresults = nsigresults+1; 
        if use_rotations==false
            fprintf('MFs outperform %s in %s\n',test_which,labels{tt})
        else
            fprintf('%s outperform their rotations in %s\n',test_which,labels{tt})
        end
    elseif (realDiff(tt)>0) && (pVals(tt)<p_corr) % mean empirical RE is significantly larger than mean control RE
        nsigresults = nsigresults+1; 
        if use_rotations==false
            fprintf('%s outperforms MFs in %s\n',test_which,labels{tt})
        else
            fprintf('%s outperforms MFs in %s\n',test_which,labels{tt})
        end
    end
end
if nsigresults==0
    fprintf('No significant difference between MFs and %s\n',test_which)
end

% all tasks together
realDiff_alltasks = mean(set2(:)) - mean(set1(:)); % hypothesis that is tested 
[pVal_alltasks, ~] =  checkDiffSignificanceMonteCarlo([set2(:),set1(:)], nperm);
if (realDiff_alltasks<0) && (pVal_alltasks<alpha) % mean empirical RE is significantly smaller than mean control RE
    if use_rotations==false
        fprintf('MFs outperform %s\n',test_which)
    else
        fprintf('%s outperform their rotations\n',test_which)
    end
elseif (realDiff_alltasks>0) && (pVal_alltasks<alpha) % mean empirical RE is significantly larger than mean control RE
    if use_rotations==false
        fprintf('%s outperforms MFs\n',test_which)
    else
        fprintf('%s outperforms MFs\n',test_which)
    end
else
    fprintf('No significant difference between MFs and %s\n',test_which)
end
%% plot correlation for each task type 
save_now = true; 
if strcmp(test_which,'rotations') && ~exist('corr_match_rand_full','var')
    corr_match_rand_full = corr_match_rand;
    corr_match_rand = min(corr_match_rand,[],3);
end
f=figure;
f.Units = 'centimeters';
for tt=1:7
    subplot(2,4,tt)
    set(gca,'FontName','FreeSans')
    hold on
    h = gobjects(2,1);
    transparent = 1;
    
    % solid curve: 
    if strcmp(test_which,'ICs15')
        this_use_bf = use_bf_emp.use_bf(matches2); 
        x = this_use_bf(2:end)-1; % non-constant harmonics
        this_task_type = corr_match(2:end,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',clrs(1,:),'LineWidth',2},transparent);
    elseif (strcmp(test_which,'harmonics') && (use_rotations==true)) || (~strcmp(test_which,'harmonics') && (use_rotations==false)) % first curve is harmonics in its color
        x = use_bf(2:end)-1; % non-constant harmonics
        this_task_type = corr_match(2:end,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',clrs(1,:),'LineWidth',2},transparent);
    elseif ~strcmp(test_which,'harmonics') && (use_rotations==false) % first curve is other basis set in its color
        x = use_bf(1:end); % first component is not constant and is used 
        this_task_type = corr_match(:,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',clrs(1,:),'LineWidth',2},transparent);
    end
    h(1) = H.mainLine;
    
    % dashed curve 
    if (strcmp(test_which,'harmonics') && (use_rotations==true)) % second curve is harmonics' rotations, in grey 
        x = use_bf(2:end)-1; % non-constant harmonics
        this_task_type = corr_match_rand(2:end,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',grayclr,'LineWidth',2,'LineStyle','--'},transparent);
    elseif ~strcmp(test_which,'harmonics') % second curve is other basis set in its color
        x = use_bf(1:end);
        this_task_type = corr_match_rand(:,tasktypes_red==tt);% & (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))');
        y = mean(this_task_type,2);
        errBar = [max(this_task_type,[],2)-y,y-min(this_task_type,[],2)];
        H = shadedErrorBar(x',y,errBar,{'Color',this_clr,'LineWidth',2,'LineStyle','--'},transparent);
    end
    
    h(2) = H.mainLine;    
%     title(labels{tt})
    if tt==5 || tt==1
        ylabel('Correlation')
    end
    if strcmp(test_which,'ICs15')
        set(gca,'XLim',[1 15])
    else
        set(gca,'XLim',[1 100])
        set(gca,'XTick',[20:40:100])
        set(gca,'XTickLabel',[20:40:100])
    end    
    set(gca,'YLim',[0 1])     
    xlabel('Components')
    set(gca,'FontSize',fs)
    grid on
end
subplot(2,4,8)
if ~strcmp(test_which,'harmonics')
    plot(1,1,'Color',clrs(1,:),'LineWidth',2);
else
    plot(1,1,'Color',this_clr,'LineWidth',2);
end
hold on
if use_rotations==true && strcmp(test_which,'harmonics')==true
    plot(1,1,'Color',grayclr,'LineWidth',2,'LineStyle','--');
else
    plot(1,1,'Color',this_clr,'LineWidth',2,'LineStyle','--');
end

if use_rotations==false % compare harmonics to controls 
    if strcmp(test_which,'eigvecs')
        [~,icons,~,~] = legend('Fct. harmonics','Adj. eigenvecs','Location','Best');
    elseif strcmp(test_which,'FC_eigvecs')
        [~,icons,~,~] = legend('Funtional harmonics','FC eigenvecs','Location','Best');
    elseif strcmp(test_which,'PCs')
        [~,icons,~,~] = legend('Funtional harmonics','PCs','Location','Best');
    elseif strcmp(test_which(1:3),'ICs')
        [~,icons,~,~] = legend('Funtional harmonics','ICs','Location','Best');
    end
elseif use_rotations==true % compare harmonics and controls to their rotations
    if strcmp(test_which,'harmonics')
        [~,icons,~,~] = legend('Funtional harmonics','Rotations','Location','Best');
    elseif strcmp(test_which,'eigvecs') && use_rotations==true
        [~,icons,~,~] = legend('Adj. eigenvecs','Rotations','Location','Best');
    elseif strcmp(test_which,'FC_eigvecs')
        [~,icons,~,~] = legend('FC eigenvecs','Rotations','Location','Best');
    elseif strcmp(test_which,'PCs')
        [~,icons,~,~] = legend('PCs','Rotations','Location','Best');
    end
end
legend boxoff
set(gca,'FontSize',fs)
% adjust legend line lengths
leglnspos = icons(3).XData;
for ch=(2+1):2:(2*3)
    icons(ch).XData = [leglnspos(2)-0.12 leglnspos(2)];
end
axis off
set(gcf,'Color','white')
f.Position = [26.3790   19.5527   0.75*twocolwidth   8];

if save_now==true && use_rotations==false
    export_fig(sprintf('SI_fig_corr_tasks_%s.png',test_which),'-r300')
elseif save_now==true && use_rotations==true
    export_fig(sprintf('SI_fig_corr_tasks_%s_rotated.png',test_which),'-r300')
end

%% plot all basis sets together
% Figure from main manuscript 
save_now = false;

zoom = [false,true]; % make one zoomed in one non-zoomed in version 
% use colors that are different from standard colors (those signify task
% groups) 
% clrs = [230,25,75;128,128,0;0,128,128;170,110,40]/255; % red, olive, teal,brown
twocolwidth = 18.3; % cm 
if exist('RE_corrs_harmonics.mat','file') && exist('RE_corrs_harmonics_rotated.mat','file') && exist('RE_corrs_eigvecs.mat','file') &&...
        exist('RE_corrs_PCs.mat','file') && exist('RE_corrs_FC_eigvecs.mat','file') && exist('RE_corrs_ICs15.mat','file')
    harmonics = load('RE_corrs_harmonics.mat');
    harmonics_rotated = load('RE_corrs_harmonics_rotated.mat'); 
    adj_eigvecs = load('RE_corrs_eigvecs.mat');
    FC_eigvecs = load('RE_corrs_FC_eigvecs.mat');
    PCs = load('RE_corrs_PCs.mat');
    ICs15 = load('RE_corrs_ICs15.mat');
    
    % plot one line per set, one subplot for each task group
    f=figure;
    f.Units = 'centimeters';
    % position/size
    pos = f.Position;
    f.Position = [pos(1) pos(2) twocolwidth/5*3.7 7.5];
    
    for tt=1:7
        % params for subtightplot
        gapvert = 0.2;
        gaphorz = 0.05;
        margin_height = [0.12 0.1]; % lower, upper
        margin_width = [0.1 0.05]; % left, right
        lw = 2;
        try
            subtightplot(2,4,tt,[gapvert gaphorz], margin_height, margin_width);
        catch
            error('Error using subtightplot. Maybe you do not have this package. Download from https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot.')
        end
        set(gca,'FontName','FreeSans')
        hold on
        % rotations
        % use best RE
        RE_match_rand = min(harmonics_rotated.RE_match_rand,[],3);
        this_task_type = RE_match_rand(2:end,tasktypes_red==tt);
        x = harmonics_rotated.use_bf;
        x = x(2:end); % leave out constant harmonic
        x = x-1; %...but plot from 1:100
        y = mean(this_task_type,2);
        h2 = plot(x,y,'Color',grayclr,'LineWidth',lw);
        % adj eigvecs
        this_task_type = adj_eigvecs.RE_match(:,tasktypes_red==tt);
        x = adj_eigvecs.use_bf;
        y = mean(this_task_type,2);
        h3 = plot(x,y,'Color',clrs(2,:),'LineWidth',lw);
        % FC eigvecs
        this_task_type = FC_eigvecs.RE_match(:,tasktypes_red==tt);
        x = FC_eigvecs.use_bf;
        y = mean(this_task_type,2);
        h4 = plot(x,y,'Color',clrs(3,:),'LineWidth',lw);
        % PCs
        this_task_type = PCs.RE_match(:,tasktypes_red==tt);
        x = PCs.use_bf;
        y = mean(this_task_type,2);
        h5 = plot(x,y,'Color',clrs(4,:),'LineWidth',lw);
        % ICs
        this_task_type = ICs15.RE_match(:,tasktypes_red==tt);
        x = ICs15.use_bf;
        y = mean(this_task_type,2);
        h6 = plot(x,y,'Color','k','LineWidth',lw);
        % harmonics
        this_task_type = harmonics.RE_match(2:end,tasktypes_red==tt);
        x = harmonics.use_bf;
        x = x(2:end); % leave out constant harmonic
        x = x-1; %...but plot from 1:100
        y = mean(this_task_type,2);
        h1 = plot(x,y,'Color',clrs(1,:),'LineWidth',lw);
        % params
        set(gca,'XLim',[0 11])%01])
        %set(gca,'XTick',[1 11])
        set(gca,'YLim',[0.5 1.5])
        xlabel('Components')
        if tt==1 || tt==5
            ylabel({'Normalized';'reconstruction error'})
        end
        t = title(labels{tt});
        pos = t.Position;
        t.Position = [pos(1) pos(2)+0.1 pos(3)];

        set(gca,'FontSize',fs)
        grid on
        
    end
    set(gcf,'Color','white')
    % legend in extra subplot
    subplot(2,4,8)
    plot(1,1,'Color',clrs(1,:),'LineWidth',2);
    hold on
    plot(1,1,'Color',grayclr,'LineWidth',2);
    plot(1,1,'Color',clrs(2,:),'LineWidth',2);
    plot(1,1,'Color',clrs(3,:),'LineWidth',2);
    plot(1,1,'Color',clrs(4,:),'LineWidth',2);
    plot(1,1,'Color','k','LineWidth',2);
    [l,icons,~,~] = legend('Functional harmonics','Fct. harm. rotations','Adjacency eigenvectors','FC eigenvectors','Principal components','Independent components','Location','Best');
    legend boxoff
    l.Position = [0.71 0.1652 0.2837 0.2308];
    set(l,'FontName','FreeSans')
    % adjust legend line lengths
    leglnspos = icons(7).XData;
    for ch=(6+1):2:length(icons)
        icons(ch).XData = [leglnspos(1)+0.15 leglnspos(2)];
    end
    for ch=1:6
        icons(ch).FontSize = 7;
    end
    axis off
    % position/size
    pos = f.Position;
    f.Position = [pos(1) pos(2) twocolwidth/5*3.7 7.5];

    if save_now==true
        export_fig('REs_all_controls.png','-r300')
    end
else
    warning('Not all files exist.')
end










