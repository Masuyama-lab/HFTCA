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
function [Model] = TrainFTCA(Samples,net,Level,NumSteps)

% TRAINFTCA  Construct the FTCA network.

Model.Samples=Samples;
DATA=Samples.'; 

numNodes = net.numNodes;       % Number of nodes
weight = net.weight;           % Node position
CountNode = net.CountNode;     % Counter for each node
adaptiveSig = net.adaptiveSig; % sigma in each node

% Parameters for Topology
edge = net.edge;
NewEdgedNode = net.NewEdgedNode;
Lambda = net.Lambda;
minCIM = net.minCIM(Level);
Model.NumSteps=NumSteps; % Total number of steps
Model.Winners=zeros(1,size(Samples,2));

for sampleNum = 1:NumSteps

    index = mod(sampleNum, size(DATA, 1));
    if index == 0
        index = size(DATA, 1); 
    end

    % Calculate an initial kernel bandwidth for CIM based on data points.
    if isempty(weight) == 1 || mod(index, Lambda) == 0
        estSig = SigmaEstimation(DATA, index, Lambda);
    end
    
    
    % Current data sample.
    input = DATA(index,:);
    
    
    %%
    if size(weight,1) < 2 % In the case of the number of nodes in the entire space is small.
        % Add Node
        numNodes = numNodes + 1;
        weight(numNodes,:) = input;
        CountNode(numNodes) = 1;
        NewEdgedNode(1, numNodes) = 0;
        edge(numNodes, :) = 0;
        edge(:, numNodes) = 0;
        adaptiveSig(numNodes) = estSig;
        
        %================================
        Model.Winners(sampleNum)=numNodes;
        
    else
        
        % Calculate CIM based on global mean adaptiveSig.
        globalCIM = CIM(input, weight, mean(adaptiveSig));

        % Set CIM state between the local winner nodes and the input for Vigilance Test.
        [Lcim_s1, s1] = min(globalCIM);
        globalCIM(s1) = inf;
        [Lcim_s2, s2] = min(globalCIM);
    
        if Lcim_s1 > minCIM
            % Add Node
            numNodes = numNodes + 1;
            weight(numNodes,:) = input;
            CountNode(numNodes) = 1;
            NewEdgedNode(1, numNodes) = 0;
            edge(numNodes, :) = 0;
            edge(:, numNodes) = 0;
            adaptiveSig(numNodes) = mean(adaptiveSig);
            Model.Winners(sampleNum)=numNodes;

        else
            weight(s1,:) = weight(s1,:) + (1/CountNode(s1)) * (input - weight(s1,:));
            CountNode(s1) = CountNode(s1) + 1;
            Model.Winners(sampleNum)=s1;

            if Lcim_s2 <= minCIM
                % Update weight of s1 neighbors.
                s1Neighbors = find(edge(s1,:));
                for k = s1Neighbors
                    weight(k,:) = weight(k,:) + (1/(10*CountNode(k))) * (weight(s1,:) - weight(k,:));
                end

                % Create an edge between s1 and s2 nodes.
                NewEdgedNode(1,s1) = 1;
                edge(s1,s2) = 1;
                edge(s2,s1) = 1;

            end
        end
    end

    
    % Topology Reconstruction
    if mod(sampleNum, Lambda) == 0 && ~isempty(weight)
        % -----------------------------------------------------------------
        % Delete Node based on number of neighbors
        nNeighbor = sum(edge);
        deleteNodeEdge = (nNeighbor == 0);
        
        % Delete process
        numNodes = numNodes - sum(deleteNodeEdge);
        weight(deleteNodeEdge, :) = [];
        CountNode(deleteNodeEdge) = [];
        NewEdgedNode(:, deleteNodeEdge) = [];
        edge(deleteNodeEdge, :) = [];
        edge(:, deleteNodeEdge) = [];
        adaptiveSig(deleteNodeEdge) = [];
        
        deleteIdx = find(deleteNodeEdge);
        for i = 1:size(deleteIdx,2)
            
            Model.Winners(Model.Winners==deleteIdx(i)) = -1;
            Model.Winners(Model.Winners>deleteIdx(i)) = Model.Winners(Model.Winners>deleteIdx(i))-1;
        end
        % -----------------------------------------------------------------
        % Delete Intersections of edge
        [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));
        
    end % if mod(sampleNum, Lambda) == 0
    
    
end % for sampleNum = 1:size(DATA,1)



% -----------------------------------------------------------------

if ~isempty(weight)
    % Delete Node based on number of neighbors
    nNeighbor = sum(edge);
    deleteNodeEdge = (nNeighbor == 0);

    % Delete process
    numNodes = numNodes - sum(deleteNodeEdge);
    weight(deleteNodeEdge, :) = [];
    CountNode(deleteNodeEdge) = [];
    NewEdgedNode(:, deleteNodeEdge) = [];
    edge(deleteNodeEdge, :) = [];
    edge(:, deleteNodeEdge) = [];
    adaptiveSig(deleteNodeEdge) = [];
    deleteIdx = find(deleteNodeEdge);
    for i = 1:size(deleteIdx,2)

        Model.Winners(Model.Winners==deleteIdx(i)) = -1;
        Model.Winners(Model.Winners>deleteIdx(i)) = Model.Winners(Model.Winners>deleteIdx(i))-1;
    end

    % -----------------------------------------------------------------
    % Delete Intersections of edge
    [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));

end

% Cluster Labeling based on edge (Functions are available above R2015b.)
connection = graph(edge ~= 0);
LebelCluster = conncomp(connection);



Model.numNodes = numNodes;      % Number of nodes
Model.weight = weight;          % Mean of nodes
Model.CountNode = CountNode;    % Counter for each node
Model.adaptiveSig = adaptiveSig;
Model.LebelCluster = LebelCluster;
Model.edge = edge;
Model.NewEdgedNode = NewEdgedNode;
Model.Lambda = Lambda;

end


% Calculate an initial kernel bandwidth for CIM based on data points.
function estSig = SigmaEstimation(DATA, sampleNum, Lambda)

if size(DATA,1) < Lambda
    exNodes = DATA;
elseif (sampleNum - Lambda) <= 0
    exNodes = DATA(1:Lambda,:);
elseif (sampleNum - Lambda) > 0
    exNodes = DATA( (sampleNum+1)-Lambda:sampleNum, :);
end

% Scaling [0,1]
normalized = (exNodes-min(exNodes))./(max(exNodes)-min(exNodes));
qStd = std(normalized);
qStd(isnan(qStd))=0;

% normal reference rule-of-thumb
% https://www.sciencedirect.com/science/article/abs/pii/S0167715212002921
[n,d] = size(exNodes);
estSig = median( ((4/(2+d))^(1/(4+d))) * qStd * n^(-1/(4+d)) );

end


% Correntropy induced Metric (Gaussian Kernel based)
function cim = CIM(X,Y,sig)
% X : 1 x n
% Y : m x n
[n, att] = size(Y);
g_Kernel = zeros(n, att);

for i = 1:att
    g_Kernel(:,i) = GaussKernel(X(i)-Y(:,i), sig);
end

ret0 = GaussKernel(0, sig);
ret1 = mean(g_Kernel, 2);

cim = sqrt(ret0 - ret1)';
end

function g_kernel = GaussKernel(sub, sig)
g_kernel = exp(-sub.^2/(2*sig^2));
% g_kernel = 1/(sqrt(2*pi)*sig) * exp(-sub.^2/(2*sig^2));
end



% Delete intersections of edge
function [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, sigma)

% for d = 1:size(weight,1); % Search all nodes
for d = find(NewEdgedNode == 1) % Search only new edged nodes
    
    node1 = find(edge(d,:)); % Neighbors of d-th node
    if size(node1,1) >= 1
       posX1 = weight(d,:); % position of d-th node
        for m = 1:size(node1,2) % Search all neighbors of d-th nodes
            posY1 = weight(node1(m),:); % position of m-th neighbor node of d-th node
            for h = 1:size(node1,2)
                target2 = node1(h);
                node2 = find(edge(target2,:)); % Neighbors of m-th node
                posX2 = weight(target2,:); % position of h-th neighbor node of m-th node
                for k = 1:size(node2,2)
                    posY2 = weight(node2(k),:); % position of k-th neighbor node of h-th node
                    isConvex = findIntersection(posX1, posY1, posX2, posY2); % find intersections
                    if isConvex == 1 % If intersection is exist, delete edge which has larger CIM.
                        cim1 = CIM(weight(d,:), weight(node1(m),:), sigma);
                        cim2 = CIM(weight(target2,:), weight(node2(k),:), sigma);
                        if cim2 >= cim1
                            edge(target2, node2(k)) = 0;
                            edge(node2(k), target2) = 0;
                        else
                            edge(d, node1(m)) = 0;
                            edge(node1(m), d) = 0;
                        end
                    end % end isConvex
                end % end k
            end % end h
        end % end m  
    end

end % end d

NewEdgedNode = zeros(size(NewEdgedNode));

end

% Check intersection of edges
function [isConvex] = findIntersection(A, B, C, D)

F1  = B(:,1)-D(:,1);
F2  = B(:,2)-D(:,2);
M11 = B(:,1)-A(:,1);
M21 = B(:,2)-A(:,2);
M12 = C(:,1)-D(:,1);
M22 = C(:,2)-D(:,2);
deter = M11.*M22 - M12.*M21;
lambda = -(F2.*M12-F1.*M22)./deter;
gamma = (F2.*M11-F1.*M21)./deter;

% E = (lambda*[1 1]).*A + ((1-lambda)*[1 1]).*B;
% isConvex = (0 <= lambda & lambda <= 1)  & (0 <= gamma & gamma <= 1);

isConvex = (0 < lambda & lambda < 1)  & (0 < gamma & gamma < 1) ;
isConvex = isConvex';

end
