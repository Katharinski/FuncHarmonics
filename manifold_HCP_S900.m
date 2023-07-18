% script to extract manifolds from HCP dense connectome of 820
% subjects-release; code adapted from Selen Atasoy by Katharina Glomb, 
% November 2016 
% katharina.glomb@gmail.com 

clearvars
global cdata % to save memory

%% parameters for manifold calculation
m_params.knn = 300; % number of nearest neighbors
m_params.dims = 100; % number of eigenvectors (+1 for constant)
method = 'nn'; % nearest neighbor search
laplacian = 'combinatorial';

%% compute the adjacency matrix, searching for nearest neighbors
% this will take several days if not parallelized - make sure to save it! 
adj_save_name = ['adjacency_nn',num2str(m_params.knn),'.mat'];
if ~exist(adj_save_name,'file')
    fprintf('Loading distance vector...\n')    
    load('HCP_S900_pdist_vec','cdata') % extracted from HCP dense FC, see extract_dconn.m
    A = computeAdjacencyFastEfficient(cdata, method, m_params.knn, 0);
    save(adj_save_name,A,m_params)
else
    load(adj_save_name,'A')
end

%% confirm that there is only one connected component
[~,no_comps] = getConnectedComponents(A);

%% compute the functional harmonics
% the functional harmonics are the columns of posCORR.M, which is an (m_params.dims+1) x (number of vertices) matrix
% first column=constant; functional harmonics shown in the paper: columns 2 to 12
if no_comps==1
    [posCORR.M, posCORR.eigVals]  = computeLaplacianEigenmaps(A, m_params.dims, laplacian);
    posCORR.comp_inds = ones(1,size(A,1));
    save(sprintf('HCP_S900_CORR_manifold_knn%i',m_params.knn), 'posCORR', 'm_params','-v7.3');
else
    warning('There was more than 1 connected component!')
end


