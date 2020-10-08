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
function ImgPlot = GetPrototypesImg(Centroids,Winners,ImgSize)

% GETPROTOTYPESIMG  Create quatization image.
% Inputs:
%   Centroids = Prototypes of the HFTCA model
%   Winners = Winning neurons for each input sample
%   ImgSize = Image size
% Output:
%   ImgProt = Quatization image

% Get the prototypes of the HFTCA model
Prototypes = Centroids(:,Winners);

% Prototypes redimension
ImgPlot = reshape(Prototypes(1:3,:), ImgSize(3), ImgSize(1), ImgSize(2));
ImgPlot = shiftdim(ImgPlot,1);