% 
% (c) 2020 Naoki Masuyama
% 
% These are the codes of Hierarchical Fast Topological CIM-based ART (HFTCA)
% proposed in "Y. Yamada, N. Masuyama, N. Amako, Y. Nojima, C. K. Loo, and H. Ishibuchi,
% Divisive Hierarchical Clustering Based on Adaptive Resonance Theory,
% Proc. of 2020 International Symposium on Community-centric Systems (CcS 2020), 
% pp. 1-6, Tokyo, Japan, September 23-26, 2020."
% 
% Please contact "masuyama@cs.osakafu-u.ac.jp" if you have any problems.
%   


clear all; close all;

% Load data
ImgOriginal = imread("baboon.tiff");

% Resize Image and reshape data
Compression_rate = 0.125;
ImgOriginal = imresize(ImgOriginal,Compression_rate);
ImgOriginal = double(ImgOriginal)/255;
Samples = reshape(shiftdim(ImgOriginal,2),3,[]);

% Randamize data
ran = randperm(size(Samples,2));
ranSamples = Samples(:,ran);

%   Parameters of HFTCA ===================================================
HFTCAnet.numNodes    = 0;                       % Number of clusters
HFTCAnet.weight      = [];                      % Mean of cluster
HFTCAnet.CountNode = [];                        % Counter for each node
HFTCAnet.edge = [];                             % Edge connections matrix
HFTCAnet.NewEdgedNode = [];                     % Node which creates new edge
HFTCAnet.adaptiveSig = [];                      % Sigma in each node
HFTCAnet.minCIM = [0.6 0.35 0.3];               % Similarity threshold
HFTCAnet.MaxLevel = size(HFTCAnet.minCIM,2);    % Max layer level
HFTCAnet.Lambda = 100;                          % Interval for node deletion and topology construction
HFTCAnet.Epochs=1;                              % Mean of cluster
%==========================================================================

% HFTCA Training
tic
[Model] = TrainHFTCA(ranSamples,HFTCAnet,1);
Time_hftca = toc;


% Plot network structure
Model2 = Generate_Model(Model); % Change the order of model nodes in hue order
[fig, num_nodes_hftca] = PlotHFTCA_GetNumNodes(Model2);
hold off

if ~isempty(Model)

    Centroids_hftca = GetCentroidsHFTCA(Model,0);
    
    % Count the number of leaf nodes
    NumLeaves_hftca = length(Centroids_hftca);

    if ~isempty(Centroids_hftca)
        Winners_hftca = GetWinnersHFTCA(Centroids_hftca, Samples);
        ImgPlot_hftca = GetPrototypesImg(Centroids_hftca, Winners_hftca, size(ImgOriginal));
        figure
        imshow(ImgPlot_hftca);
        hold off

        % Compute the PSNR
        psnr_r = PSNR(ImgOriginal(:,:,1), ImgPlot_hftca(:,:,1));
        psnr_g = PSNR(ImgOriginal(:,:,2), ImgPlot_hftca(:,:,2));
        psnr_b = PSNR(ImgOriginal(:,:,3), ImgPlot_hftca(:,:,3));
        PSNR_hftca = (psnr_r+psnr_g+psnr_b)/3;

        % Compute the Davies-Bouldin Index
        DBI_hftca = db_index(Samples', Winners_hftca, Centroids_hftca');

        % Disp results
        fprintf('\nResults ================\n')
        fprintf('Time_hftca: %d\n',Time_hftca);
        fprintf('PSNR_hftca: %d\n',PSNR_hftca);
        fprintf('DBI_hftca: %d\n',DBI_hftca);
        fprintf('NumLeaves_hftca: %d\n',NumLeaves_hftca);
    end
    
end 
