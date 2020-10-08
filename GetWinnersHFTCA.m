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
function [Winners, Errors] = GetWinnersHFTCA(Prototypes,Samples)

% Determine the winner node based on distance
%
% Inputs:
%   Prototypes = Prototypes of the HFTCA model
%   Samples = Matrix of input data
% Outputs:
%   Winners = Winning neurons for each input sample
%   Errors = Errors for each input sample

NumSamples = size(Samples,2);
Winners = zeros(NumSamples,1);
Errors = zeros(NumSamples,1);
NumNeuro = size(Prototypes,2);

% Main loop
for SampleIndex = 1:NumSamples
    SquaredDistances = sum((repmat(Samples(:,SampleIndex),1,NumNeuro)-Prototypes).^2,1);
    [MinimumError,WinnerIndex] = min(SquaredDistances);
    Winners(SampleIndex) = WinnerIndex;
    Errors(SampleIndex) = MinimumError;    
end