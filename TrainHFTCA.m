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
function [Model] = TrainHFTCA(Samples,net,Level)

% TRAINFTCA  Create the FTCA tree.

%%
Model = [];
MaxLevel = net.MaxLevel;
[Dimension,NumSamples]=size(Samples);
if ((NumSamples<(Dimension+1)) && (Level>1)) || (Level>MaxLevel)    
    return;
end

%fprintf('\nLEVEL=%d\n',Level);

%% Growing Process
NumSteps = net.Epochs*NumSamples;
Model = TrainFTCA(Samples,net,Level,NumSteps);

%% Expansion Process
Winners = Model.Winners;
Model.Means = Model.weight';
NeuronsIndex = find(isfinite(Model.Means(1,:)));
NumNeurons = numel(NeuronsIndex);

%fprintf('Final Graph Neurons: %d\n',NumNeurons);

%%
Model.Connections = sparse(Model.edge);

%% PRUNE THE GRAPHS WITH ONLY 2 NEURONS. THIS IS TO SIMPLIFY THE HIERARCHY
if NumNeurons==2
    Model=[];
    return;
else
    for NeuronIndex=NeuronsIndex
        ChildSamples = Samples(:,Winners==NeuronIndex);
        Model.Child{NeuronIndex} = TrainHFTCA(ChildSamples,net,Level+1);
    end
end

end

