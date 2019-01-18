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


%% Run Graph Clustering with Symmetric NMF Algorithm
cN = 4; % We choose 4 cluster number (see paper for how this was chosen, 
%code for the estimation is shown below)

clus_weights = symnmf_anls(simMatrix,cN);
% clus_weights is the cluster weights of each activation map
% We chose 4 as the cluster number 

%% Visualize Cluster Centroids in Nifti Space 

% In order to replicate the centroid figures, run the code presented in
% this section. We first take a weighted average centroid (based off cluster
% weights). Next we project this weighted centroid onto voxel space using
% the H weights of the NMF solution computed earlier. Because of the
% sub-sampling of the activation maps before computation of the NMF
% (described in the paper), we have to use an interpolation procedure to
% produce a full nifti map for visualization.
    
    % 1. Compute Centroid and Project to Voxel Space
    centroid_voxel_weights = zeros(size(H,2),com);
    centroid_weights = zeros(size(W,2),com);
    for i = 1:com
       weights = clus_weights(:,i); 
       for j = 1:size(W,2)
           temp = W(:,j);
           centroid_weights(j,i) = wmean(temp,weights);
       end
       centroid_voxel_weights(:,i) = centroid_weights(:,i)'*H;
    end 
    
    % 2. Interpolate in Nifti Space
       % i) Subsample Standard MNI Template to get Indices
           %First, we have to get the subsampled indices of the nifti space
           % (this was used previously to subsample the original activation
           % maps).
           TEMPLATE = load_nii('/Users/tbolt/Desktop/University of Miami/BCCL Research/BrainMap Project/MACM15/Grey10.nii');

            %%%% Load reduced data
            xleer_reduced      = TEMPLATE.img>0.1;
            xleer2_reduced     = single(zeros(size(TEMPLATE.img)));

            xleer2_reduced(1:2:end,1:2:end,1:2:end) = 1;
            GMindices_reduced = find(xleer_reduced.*xleer2_reduced);

            mat = vertcat(TEMPLATE.hdr.hist.srow_x,TEMPLATE.hdr.hist.srow_y...
                ,TEMPLATE.hdr.hist.srow_z);

            [X Y Z] = ind2sub(size(TEMPLATE.img),find(xleer_reduced.*xleer2_reduced));
            GMxyz_reduced = [X Y Z]'; clear X Y Z
            GMmm_reduced  = mat * [GMxyz_reduced; ones(1,size(GMxyz_reduced,2))];
            GMmm_reduced  = GMmm_reduced(1:3,:)';
            clear xleer_reduced xleer2_reduced

            %%%% Load full data
            xleer_full = TEMPLATE.img>0.1;
            xleer2_full = single(zeros(size(TEMPLATE.img)));

            xleer2_full(1:end,1:end,1:end) = 1;
            GMindices_full = find(xleer_full.*xleer2_full);

            [X Y Z] = ind2sub(size(TEMPLATE.img),find(xleer_full.*xleer2_full));
            GMxyz_full = [X Y Z]'; clear X Y Z
            GMmm_full  = mat * [GMxyz_full; ones(1,size(GMxyz_full,2))];
            GMmm_full  = GMmm_full(1:3,:)';
            clear xleer_reduced xleer2_reduced
            
        % ii) Put into 3d Nifti Space and Interpolate using 'inpaintn'
        % function
            centroid_maps = zeros(91,109,91,size(W,2));
            for i=1:size(centroid_voxel_weights,2)
                  disp(i)
                  map_centroid = squeeze(centroid_maps(:,:,:,i));
                  map_centroid(GMindices_full) = NaN;
                  map_centroid(GMindices_reduced) = centroid_voxel_weights(:,i);
                  interp_map_centroid = inpaintn(map_centroid);
                  centroid_maps(:,:,:,i) = interp_map_centroid;
            end
            
    % 3. Create Niftis using the 'Nifti' Toolbox
    nii = make_nii(centroid_maps,[2 2 2]);
    save_nii(nii,strcat('Centroid_Activation_Clusters_SymNMF_2.nii'))

%% Multivariate Distance Matrix Regression Analysis

    % 1. Dummy-Code BrainMap Descriptors
    % Because of data sharing issues we only supply the already dummy-coded
    % matrices. The 

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

%% Optional Code: Split-Half Estimation of the Number of Clusters for the Symmetric NMF Algorithm
% I would use parallel processing toolbox to speed up the inner loop in this process
% Pre-allocate the cluster 'instability' measure, which we use to 
% choose the number of clusters
instability = zeros(20,100);

% We loop through a cluster solution of 1 to 20
for i = 1:20
    disp(i)
    for j = 1:100
        idx = randperm(numel(1:8919)); % Get a random permutation of indices to switch around observations
        sample1 = double(W(idx(1:4459),:)); % Split the data into two random halves
        sample2 = double(W(idx(4460:end),:));
        
        simMatrix1 = 1 - squareform(pdist(sample1,'cosine')); % Compute Similarity Matrix for split-half
        H1 = symnmf_anls(simMatrix1,i); % Run Sym-NMF algorithm on split-half 
        simMatrix2 = 1 - squareform(pdist(sample2,'cosine'));
        H2 = symnmf_anls(simMatrix2,i);
        
        % Compute the centroids of each split-half using a weighted average
        centroid_weights1 = zeros(size(W,2),i);
        centroid_weights2 = zeros(size(W,2),i);
        for c = 1:i
           weights1 = H1(:,c);
           weights2 = H2(:,c);
           for h = 1:size(W,2)
               temp1 = sample1(:,h);
               temp2 = sample2(:,h);
               centroid_weights1(h,c) = sum(temp1.*weights1)./sum(weights1);
               centroid_weights2(h,c) = sum(temp2.*weights2)./sum(weights2);
           end
        end 
        % Use Amari-Max Error to measure 'instability' between the split-half cluster
        % solutions based off the correlations between their centroids
        CORR = corr(centroid_weights1,centroid_weights2);                
        instability(i,j) = amariMaxError(CORR);
        

    end
end
