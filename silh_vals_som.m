% compute somatotopy index for testing correspondence between
% functional harmonics (+ control basis function sets) and somatotopic sub-
% areas defined in Glasser et al. 2016 
%
% Katharina Glomb 
% last modified June 24, 2021
% katharina.glomb@gmail.com 

clearvars

% compute silhouette values 
test_which = 'harmonics'; % harmonics, eigvecs, PCs
suffix = ''; % optional suffix for testing purposes 
path_to_surfaces = [pwd,'/data'];

if ~exist('HCP_plot_labels_w_somatotopy.mat','file')
    error('File with surface labels of HCP parcellation including somatototopic areas does not exist. Get it by following instructions in create_som_subareas.m')
else
    load('HCP_plot_labels_w_somatotopy','surface_labels')
end

%% set up file names 
switch test_which 
    case 'harmonics'
        load('HCP_S900_CORR_manifold_knn300.mat') 
%         use_dims = (2:12); % discard constant harmonic
        use_dims = [4,8,12]; % only look at dimensions that have obvious som features
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
end

save_fname_emp = ['silh_coeffs_',test_which,'_som.mat'];
save_fname_rand = ['silh_coeffs_spherical_',test_which,'_som.mat'];

%% compute non-randomized silh vals

if ~exist(save_fname_emp,'file')
    fprintf('Computing empirical silhouette coefficients...\n')
    % for each dimension separately
    [SC_vs_rest,SC_vs_som,M_within,M_between_rest,M_between_som] = compute_silh_vals_som(posCORR.M(:,use_dims), surface_labels);
    
    % all dimensions together
%     Silh_vals_all = zeros(length(som_ROIs),1); 
%     M_within_all = zeros(length(som_ROIs),1); 
%     M_between_all = zeros(length(som_ROIs),1);
%     for n=1:numel(som_ROIs)
%         [Silh_vals_all(n),M_within_all(n),M_between_all(n)] = compute_silh_vals2(posCORR.M(:,use_dims),som_ROIs(n),neighboring_ROIs{n},surface_labels,1);
%     end
    save(save_fname_emp,'use_dims','SC_vs_rest','SC_vs_som','M_within','M_between_rest','M_between_som')
    fprintf('Done!\n')
else
    load(save_fname_emp,'SC_vs_rest','SC_vs_som','M_within','M_between_rest','M_between_som')
end

%% compute silh vals of rotations

load('Ind_S900','indices')
nV = length(indices);

alpha = 0.05; % significance level; one-sided test (silh vals assumed to be bigger in real data) 

if ~exist(save_fname_rand,'file')
    % determine number of randomizations
    nrand = (1/alpha) * length(use_dims) * 5; % 5 som regions (hems will be averaged)
    % inititalize arrays
    SC_vs_rest_rand = zeros(10,length(use_dims),nrand); % 10 bc 5 som regions on each side
    SC_vs_som_rand = zeros(10,length(use_dims),nrand);
    M_within_rand = zeros(10,length(use_dims),nrand);
    M_between_rest_rand = zeros(10,length(use_dims),nrand);
    M_between_som_rand = zeros(10,length(use_dims),nrand);
    
    % load or create rotations
    if ~exist('rotations_harmonics.mat','file')
        rotations = create_rotations(nrand, './);
        save('rotations_harmonics.mat','rotations')
    else
        load('rotations_harmonics','rotations')
        if size(rotations,2)<nrand % #rotations necessary for som. index = 5*#rotations necessary for modified silhouette value
            warning('Number of pre-computed rotations does not match indicated number of randomizations. Computing additional rotations...')
            additional_rotations = create_rotations(nrand-size(rotations,2));
            all_rotations = zeros(size(rotations,1),nrand);
            all_rotations(:,1:size(rotations,2)) = rotations;
            all_rotations(:,size(rotations,2)+1:end) = additional_rotations;
            rotations = all_rotations; 
            save('rotations_harmonics.mat','rotations')
        end
    end    
    
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
        [SC_vs_rest_rand(:,:,n), SC_vs_som_rand(:,:,n),M_within_rand(:,:,n),...
            M_between_rest_rand(:,:,n),M_between_som_rand(:,:,n)] = ...
            compute_silh_vals_som(M_rand, surface_labels);
    end
    save(save_fname_rand,'use_dims','SC_vs_rest_rand','SC_vs_som_rand','M_within_rand','M_between_rest_rand','M_between_som_rand')
    fprintf('Done!\n')
else
    load(save_fname_rand)
    nrand = size(SC_vs_rest_rand,3);
end

%% compute "somatotopy index"
% EMP
som_index_allreg = (M_between_som(:,use_dims-1)+M_between_rest(:,use_dims-1))./...
    max(M_between_som(:,use_dims-1),M_between_rest(:,use_dims-1)).*M_between_som(:,use_dims-1);
% average over right and left side
som_index = mean(cat(3,som_index_allreg(1:5,:),som_index_allreg(6:end,:)),3);

% RAND
nrand_red = nrand; 
som_index_allreg_rand = (M_between_som_rand(:,use_dims-1,1:nrand_red)+M_between_rest_rand(:,use_dims-1,1:nrand_red))./...
    max(M_between_som_rand(:,use_dims-1,1:nrand_red),M_between_rest_rand(:,use_dims-1,1:nrand_red)).*M_between_som_rand(:,use_dims-1,1:nrand_red);
% average over right and left side
som_index_rand = mean(cat(4,som_index_allreg_rand(1:5,:,:),som_index_allreg_rand(6:end,:,:)),4);

%% ...and plot
labs = {'Face','Eye','Hand','Trunk','Foot'};
labs_short = {'f','e','h','t','f'};
meanvals_dims = zeros(length(use_dims),1); 
meanvals_dims_rand = zeros(length(use_dims),nrand); 
gapvert = 0.1;
gaphorz = 0.1;
margin_height = [0.25 0.1]; % lower, upper
margin_width = [0.1 0.06]; % left, right
f=figure;
f.Units='centimeters';
pos = f.Position;
twocolwidth = 18.3; % width of two col figures in cm
f.Position = [pos(1) pos(2) twocolwidth 5.13];
clrs = get(gca,'ColorOrder');
for d=1:length(use_dims)
    subtightplot(1,length(use_dims),d,[gapvert gaphorz], margin_height, margin_width);
    %subplot(3,4,d)
    h2 = plot(squeeze(som_index_rand(:,d,:)),'x','Color',[0.5 0.5 0.5]);
    hold on
    h1 = plot(som_index(:,d),'o','LineWidth',2,'Color',clrs(1,:));
    set(gca,'XTick',1:5)
    set(gca,'FontSize',8)
    %set(gca,'YLim',[0 0.02])
    if d==1        
        ylabel('Somatotopy index')
    end
    if d==1
        legend([h1,h2(1)],{'Fct. harmonics','Rotated'},'Location','Best') 
    end
    set(gca,'XTickLabel',labs)
    RotateXLabel(90,labs,0);%0.0016);
    title(sprintf('Functional harmonic %i',use_dims(d)-1))
    grid on
    box off
end
set(gcf,'Color','white')



