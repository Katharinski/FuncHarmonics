% script for obtaining reconstructions of task maps using functional
% harmonics 
% this script 
% -computes the reconstruction errors of the harmonics - for recon. errors
% of all other basis sets and of rotations, see controls.m 
% -computes these reconstruction errors using fct. harmonics in order of
% their power in the graph spectrum
% -produces plots shown in figure 6: task maps, their reconstructions, and
% the graph Fourier spectra (spider plots), as well as confusion matrices 
% code by Katharina Glomb 
% last modified June 26, 2021 
% katharina.glomb@gmail.com 

clearvars

% load the functional harmonics 
load('HCP_S900_CORR_manifold_knn300.mat','posCORR') % posCORR.M contains fct harmonics
load('Ind_S900','indices') % indices of the medial wall 

% load the task maps 
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

neg_maps_IDs = [12,13,23:26,34,35,50:56,67,68,72,73,79,80,84,85]; % some maps are just inversions of others [?!] - remove them
id_contrasts_IDs = [14,27:30,36,57:62,66,74,78,86];
reduced_flag = false;
if reduced_flag==false
    nT_full = size(taskMaps,2);
    reduced_flag = true;
    taskMaps = taskMaps(:,(~ismember(1:nT_full,[neg_maps_IDs,id_contrasts_IDs])));
    nT = size(taskMaps,2);
end

% load cortical surfaces 
path_to_surfaces = [pwd,'/data'];
subject_name = 'S900';
% surf_type = 'inflated_MSMAll';
surf_type = 'flat';
try
    [vertices, faces] = connRSMreadGII(path_to_surfaces, subject_name, surf_type); % get surfaces
catch
    error('Error reading surface. Perhaps gifti() function is not working. Try using gifti toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
end

%% compute reconstruction errors between fct harmonics and task maps 
recon_errs_filename = 'RE_corrs_harmonics.mat';
recons_filename = 'task_recons.mat';

if ~exist(recon_errs_filename,'file')
    use_bf = [1:21,31:10:101];
    [spec,recons,RE_all,RE_match,corr_all,corr_match,p_corr_match] = compute_recon_errors(posCORR.M,taskMaps,use_bf);
    save(recon_errs_filename,'spec','RE_all','RE_match','corr_all','corr_match','p_corr_match','use_bf')
    save(recons_filename,'recons','use_bf')
else
    load(recon_errs_filename)
end

%% compare to random (rotations) 
% file needs to be created using controls.m to guarantee consistency. 
recon_errs_rand_filename = 'RE_corrs_harmonics_rotated.mat';

if ~exist(recon_errs_rand_filename,'file')
    error('File %s does not exist. Use controls.m to create it.',recon_errs_rand_filename)
else
    use_bf_rand = load(recon_errs_rand_filename,'use_bf'); % first test if correct dimensions were used 
    if ~isequal(use_bf_rand.use_bf,use_bf)
        error('use_bf not consistent between random and empirical data. Aborting.')
    else
        load(recon_errs_rand_filename)
        nrand = size(spec_rand,3);
        if nrand<((1/0.05) * 11) % we usually discuss 11 harmonics 
            warning('Number of randomizations is too low! Press any key to continue anyway.')
            pause
        end
    end
end

%% do reconstructions using the functional harmonics in the order in which they contribute 

recon_errs_sort_filename = 'RE_corrs_sort.mat';
recons_sort_filename = 'task_recons_sort.mat';

if ~exist(recon_errs_sort_filename,'file')
    use_bf = [1:21,31:10:101];
    [~,recons_sort,RE_all_sort,RE_match_sort,corr_all_sort,corr_match_sort,~] = compute_recon_errors(posCORR.M,taskMaps,use_bf,true);
    save(recon_errs_sort_filename,'spec','RE_all_sort','RE_match_sort','corr_all_sort','corr_match_sort','use_bf')
    save(recons_sort_filename,'recons_sort','use_bf')
else
    load(recon_errs_sort_filename)
end

%% plot task maps
load(recons_sort_filename,'recons_sort')
V = size(posCORR.M,1);
set(0,'DefaultAxesFontName','FreeSans')
save_now = false;
% load(recons_filename,'recons')
for t=34%1:nT
    %close all
    nrmlz = true; % scale between -1 and 1
    colconst = false;
    % plot task
    mini = min(min([taskMaps(1:V,t),recons_sort(:,:,t)]));
    maxi = max(max([taskMaps(1:V,t),recons_sort(:,:,t)]));
    fname_tasks = sprintf('taskmap%i_small',t);
    plot_HCP_on_surface(taskMaps(1:V,t),vertices,faces,nrmlz,mini,maxi,colconst);
    
    % plot reconstructions starting with highest contribution
    % save_now = true;
    this_spec = abs(spec(:,t))/max(abs(spec(:,t)));
    [~,maxID] = max(this_spec);
    if maxID==1
        this_use_bf = [2,4,23]; % if strongest contributor is constant harmonic, show second
    else
        this_use_bf = [2,5,6,10,29];%[1,4,23];
    end
    for bf=5%this_use_bf%[2:12,29]
        plot_HCP_on_surface(recons_sort(:,bf,t),vertices,faces,nrmlz,mini,maxi,colconst);
        set(gca,'FontName','FreeSans')
        if save_now
            if bf<22 && ~(maxID==1 && bf==2)
                export_fig(sprintf('taskmap%i_recon_sort_%i_small.png',t,use_bf(bf)),'-r300'); % due to weird bug, pdf can't do the font in this case
            else
                export_fig(sprintf('taskmap%i_recon_sort_%i_small.png',t,use_bf(bf)-1),'-r300'); % due to weird bug, pdf can't do the font in this case
            end
        end
    end
end

%% spider plots (graph Fourier spectra) 
save_now = false;
% get original task map ID back for saving 
all_IDs = 1:nT_full; 
IDs_red = all_IDs((~ismember(all_IDs,neg_maps_IDs)) & (~ismember(all_IDs,id_contrasts_IDs)));
ffs = 7;
twocolwidth = 18.3;
order = 2:12;%[1,2,4,3,7,11,5,6,8,9,10];
clrs = get(gca,'ColorOrder');
for t=1:nT% [1,2,11,15]
    this_spec = abs(spec(:,t))/max(abs(spec(:,t))); % normalize to largest contributor (out of 101)
    this_spec = this_spec(order); % only plot selected part of spec
    f=figure;
    N = 11;
    theta = linspace(0.,2*pi-2*pi/N,N);
    tt = tasktypes_red(t); 
    polarplot([theta,0],[abs(this_spec(1:N));this_spec(1)],'LineWidth',2,'Color',clrs(tt,:))
    set(gca,'ThetaDir','clockwise')
    set(gca,'ThetaTick',linspace(0,360-360/N,N))
    set(gca,'ThetaTickLabel',order-1)
    rlim([0,1])
    set(gca,'RTick',0:0.5:1)
    set(gca,'RTickLabel','')
    set(gcf,'color','white')
    set(gca,'FontSize',ffs)
    set(gca,'FontName','FreeSans')
    f.Units = 'centimeters';
    f.Position = [11.4375    5.8958    twocolwidth/8    twocolwidth/8];
    original_task_no = IDs_red(t); 
    if save_now
        export_fig(sprintf('spider%i_small.pdf',original_task_no),'-r300');
    end
    close all 
end

%% confusion matrices with sorted contributions

conf_matrices_sort = zeros(length(use_bf),nT,nT);
for b=1:length(use_bf)
    [maxi,maxi_ID] = min(squeeze(RE_all_sort(b,:,:)),[],2);
    %[maxi,maxi_ID] = max(squeeze(corr_all_sort(b,:,:)),[],2); % fit of reconstruction from task to all task maps is in ROWS of corr_all(b,:,:)
    for t=1:nT
        conf_matrices_sort(b,t,maxi_ID(t)) = 1;
    end
end

%% plot selected confusion matrices (main ms)
save_now = false;
fs = 8;
f = figure;
counter = 1;
for b=[find(use_bf==1),find(use_bf==4),find(use_bf==41)]%:length(use_bf)%(1:end-1))
    subplot(1,3,counter)
    set(gca,'FontName','FreeSans')
    %imagesc(squeeze(conf_matrices_sort(b, (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs])), (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs])))))
    imagesc(squeeze(conf_matrices_sort(b,:,:)))
    %title(sprintf('%i harmonics',use_bf(b)))
    colormap(flip(colormap(cold)))
    axis square
    if counter==1
        ylabel('Reconstructions')
    end
    xlabel('Task maps')
    counter = counter+1;
    if use_bf(b)==1
        title('strongest')  
    elseif use_bf(b)<22
        title(sprintf('strongest %i',use_bf(b)))
    else
        title(sprintf('strongest %i',use_bf(b)-1))
    end
    set(gca,'FontSize',fs)
    %set(gca,'XTick',[1,19,20,22,23,35,36,38,39,41,42,44,45,47])
    %set(gca,'YTick',[1,19,20,22,23,35,36,38,39,41,42,44,45,47])
    %grid on
end
f.Units = 'centimeters';
f.Position = [11.4375    5.8958    twocolwidth/5*3    twocolwidth/4.5];
set(gcf,'Color','white')
if save_now
    export_fig('confusion_matrices_sort.png','-r300');
end

%% plot number of matched maps 
save_now = false;
nmatched_sort = zeros(length(use_bf),1);
task_matched_sort = zeros(nT,length(use_bf));
for b=1:length(use_bf)
    %nmatched_sort(b) = sum(diag(squeeze(conf_matrices_sort(b, (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs])), (~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]))))));
    nmatched_sort(b) = sum(diag(squeeze(conf_matrices_sort(b,:,:))));
    for t1=1:nT
        if conf_matrices_sort(b,t1,t1)==1
            task_matched_sort(t1,b) = 1;
        else
            mtch = find(conf_matrices_sort(b,t1,:)); % matching task
            if tasktypes(t1)==tasktypes(mtch)
                task_matched_sort(t1,b) = 1;
            end
        end
    end
end
f=figure;
set(gca,'FontName','FreeSans')
h1 = plot(use_bf,nmatched_sort/(nT),'k','LineWidth',2);
hold on
%plot([12,12],[0,1],'r')
h2 = plot(use_bf,sum(task_matched_sort)/(nT),'k-');
[~,icons,~,~] = legend([h1,h2],{'Exact task','Task group'},'Location','Best','FontSize',fs);
legend boxoff
xlabel('# of fct. harmonics')
ylabel({'Proportion of';'tasks identified'})
set(gcf,'Color','white')
set(gca,'FontSize',fs)
set(gca,'XLim',[1 101])
%set(gca,'XLim',[-1 41])
set(gca,'XTick',[50 100])
%set(gca,'YLim',[min([sum(task_matched_sort(~ismember(1:nT,[neg_maps_IDs,id_contrasts_IDs]),:))/(nT-length([neg_maps_IDs,id_contrasts_IDs])),nmatched_sort'/(nT-length([neg_maps_IDs,id_contrasts_IDs]))]) 1])
% set(gca,'YLim',[min([sum(task_matched_sort)/(nT),nmatched_sort'/(nT)]) 1])
set(gca,'YLim',[0 1])
box off
f.Units = 'centimeters';
f.Position = [11.4375    5.8958    3.9952    3.4396];

% adjust legend line lengths
leglnspos = icons(3).XData;
for ch=3:(length(icons))
    icons(ch).XData = [leglnspos(2)-0.1 leglnspos(2)];
end

if save_now
    export_fig('tasks_matched_sort.png','-r300');
end