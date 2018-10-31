%% Load Data Matrix
% Load in your activation map data matrix, the dimensions should be an
% activation map (# of activation maps in analysis) * voxel (# voxels in
% your activation maps) matrix. There is a variety of ways to do this but
% it's important that the voxel dimension is consistent across activation
% maps, meaning voxel 1 (row 1) should be the same voxel across all
% activation maps
% 
% Of note, you can run this pipeline with modeled activation maps computed
% from the BrainMap database, or you can run it with any sort of
% (unthresholded or thresholded) activation maps derived from a GLM. 

load('data_matrix.mat'); % We refer to the variable loaded in as 'data_matrix' throughout
N_ActMaps = size(data_matrix,1);
N_voxels = size(data_matrix,2);
%% Choose number of dimensions/networks (K) output from the NMF algorithm
K = 20;

%% SVD Initialization for NMF Algorithm
[W,H] = NNDSVD(double(data_matrix),K,0); 

%% Run NMF with SSVD Initialization
[W,H] = nmf(data_matrix,K,'TYPE','sparse','W_INIT',W,'H_INIT',H);

% W is the network expression matrix
% H is the voxel weight matrix for each network
%% Compute Similarity Matrix from Network Expression Estimates (W)
    sim_matrix = 1 - squareform(pdist(W,'cosine')); 
% I subtract the matrix from 1 to get a similarity matrix because cosine distance varies from 0 to 1
% on POSITIVE data (modeled activation maps from BrainMap are all
% positive), if you have activation maps with positive and negative values
% this will produce data with negative and positive similarity values,
% which may be OK with you, but be aware. 

% On a side note, you can use any distance metric you want here. Euclidean
% distance is a common metric, but is affected by mean differences (w/o some
% form of normalization beforehand) in activation across the brain which 
% may not be what you want. Distances such as cosine and correlation
% account for mean differences which may be more desirable. And of course
% one can compute distances/similarity using some Kernel function (gaussian, etc.).

%% Choose Resolution Parameter for Graph Clustering
% The Louvain graph clustering algorithm from the BCT toolbox comes equipped with a
% resolution (gamma) term to control the resolution of the clustering
% (default = 1).

gamma = 1;

%% Run Graph Clustering with Louvain

[M, Q1] = community_louvain(sim_matrix,gamma);

% M is the cluster assignments of each activation map
% Q1 is the modularity value of the cluster solution


%% Visualize Cluster Solution
[X,Y,INDSORT] = grid_communities(M);
imagesc(sim_matrix(INDSORT,INDSORT));           
hold on;                                 
plot(X,Y,'r','linewidth',2); 
title('Cluster Solution')
colormap jet

%% Optional Code: Split-Half Estimation of # dimensions to choose in NMF algorithm
Min_K = 20; % Minimum number of dimensions to estimate
Max_K = 75; % Maximum number of dimensions to estimate
N_split_half = 50; % Number of split-half samples to generate
for i = Min_K:Max_K
    disp(['NMF: ' num2str(i)])
    instability = zeros(1,N_split_half);
    for j = 1:N_split_half
        idx = randperm(numel(1:N_ActMaps)); % randomize the indices of each activation map 
        sample1 = data_matrix(idx(1:floor(N_ActMaps/2)),:); % Get one side of split-half
        sample2 = data_matrix(idx((floor(N_ActMaps/2)+1):end),:); % Get other side
        
        [W,H] = NNDSVD(sample1,i,0);
        [W1,H1,iter] = nmf(sample1,i,'TYPE','sparse','W_INIT',W,'H_INIT',H);
        
        [W,H] = NNDSVD(sample2,i,0);
        [W2,H2,iter] = nmf(sample2,i,'TYPE','sparse','W_INIT',W,'H_INIT',H);
        CORR = corr(H1',H2');                
        instability(j) = amariMaxError(CORR); % calculate instability between the two split-half samples
    end
    avg_error(i) = mean(instability); % compute the average instability across split-half sample pairs
    var_error(i) = std(instability); % compute the standard deviation of instability across split-half sample pairs
end

%% Optional Code: Consensus Clustering to Find Optimal Resolution for Graph Clustering
gamma_min = 0.7; % Minimum resolution
gamma_max = 3; % Maximum resolution
step_size = 0.05; % Step size from minimum to maximum resolution
N_samples = 50; % Number of times to run Louvain algorithm to assess consensus 
gamma = [gamma_min:step_size:gamma_max];  

% Assess consistency at each value of gamma using the adjusted rand_index
Partition_Matrix = zeros(N_ActMaps,N_samples);
for i = 1:length(gamma)

    for j = 1:N_samples
        disp(strcat(num2str(gamma(i)),': ',num2str(j)))
        [M, Q1] = community_louvain(sim_matrix,gamma(i));
        Partition_Matrix(:,j) = M;
    end 

    Rand_coeff = zeros(size(Partition_Matrix,2),size(Partition_Matrix,2));
    for j = 1:size(Partition_Matrix,2)
        for h = 1:j
            [~,Rand_coeff(j,h)] = zrand(Partition_Matrix(:,j),Partition_Matrix(:,h));
        end
    end
    Rand_coeff(logical(eye(size(Rand_coeff)))) = 0;
    Rand_coeff = squareform(Rand_coeff,'tovector');
    avg_Rand(i) = mean(Rand_coeff); % Average rand index measuring the 
    %consistency across all runs of Louvain algorithm
    


end 

% Based off examination of the consistency values from the rand index,
% choose a gamma value and compute the Partition Matrix above and compute
% the consensus matrix across all runs to get your consensus clustering. 

    D = agreement(Partition_Matrix,150); % Compute agreement matrix from all runs of Louvain, see function for second input option
    D = D/size(Partition_Matrix,2); % Normalize to be between 0 to 1
    consensus(i) = consensus_und(D,0.5,50); % there are two optional parameters here; see function for definition.
