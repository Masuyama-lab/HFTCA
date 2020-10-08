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
function [NewCentroids] = GetCentroidsHFTCA(Model,A)

% GETCENTROIDSHFTCA  Get recursively the centroids of the HFTCA model.
%
% Inputs:
%   Model = HFTCA model
% Output:
%   NewCentroids = Prototypes of the HFTCA model

NewCentroids = [];

if Model.numNodes == 0
    NewCentroids = [NewCentroids A];
else
    IndexValidNeurons = find(isfinite(Model.Means(1,:)));

    for NeuronIndex=IndexValidNeurons
        if ~isempty(Model.Child{NeuronIndex})
            M = Model.Means(:,NeuronIndex);
            ChildCentroids = GetCentroidsHFTCA(Model.Child{NeuronIndex},M); 
            NewCentroids = [NewCentroids ChildCentroids];
        else
            NewCentroids = [NewCentroids Model.Means(:,NeuronIndex)];
        end    
    end
end
