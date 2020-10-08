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
function [fig,num_nodes] = PlotHFTCA_GetNumNodes(Model)

% PLOTHFTCA_GETNUMNODES  Plot HFTCA Network and Calculate the total number
% of nodes in the network.

k = 1;
Model.GraphIndex = k;
queue{1} = Model;
label(k) = size(Model.Samples,2);

Node_value(1,:) = [0 0 0]; % Set loot node value

while (~isempty(queue))
   currentModel = queue{1};
   if (length(queue) > 1)
       queue(1) = []; %pop neuron
   else
       queue = []; %pop neuron
   end

   NeuronsIndex = find(isfinite(currentModel.Means(1,:)));
   fatherIndex = currentModel.GraphIndex;

    for NeuronIndex = NeuronsIndex
        childIndex = k + 1;
         if ~isempty(currentModel.Child{NeuronIndex})
             currentModel.Child{NeuronIndex}.GraphIndex = childIndex;
             label(childIndex) = size(currentModel.Child{NeuronIndex}.Samples,2); 
             if (~isempty(queue))
                queue = [queue, currentModel.Child(NeuronIndex)];
             else
                queue = [currentModel.Child(NeuronIndex)];
             end
         else
             label(childIndex) = sum(currentModel.Winners == NeuronIndex);
         end

        s(k) = fatherIndex;
        t(k) = childIndex;
        k = k + 1;
        Node_value(childIndex,:) = currentModel.Means(:,NeuronIndex).'; % Set node value
    end

end

s = s(~isnan(s));    
t = t(~isnan(t));

% Get the number of nodes
[~, num_nodes] = size(t);


%% Set edge value
[A,B] = size(Node_value);
for i = 1:A-1
    Edge_value(i,:) = Node_value(i+1,:);
end

%% Set node size
%nodesize settings
start_size = 0.01;
lay1_size = 4;
lay2_size = 2;
lay3_size = 1.5;

%initialization
node_count = 1;
node_size(node_count) = start_size;

%lay1
lay1 = find(s==1);
lay1 = lay1+1;
[a1,b1] = size(lay1);
for i = 1:b1
    node_count = node_count+1;
    node_size(node_count) = lay1_size;
end
%lay2
lay2_start = node_count;
for j = lay2_start:k-1
    check = 0;
    for i = lay1
        if s(j) == i
            check = 1;
        end
    end
    if check == 1
        node_count = node_count+1;
        node_size(node_count) = lay2_size;
    else
        break;
    end
end
%lay3
lay3_start = node_count;
for i = lay3_start:k-1
    node_count = node_count+1;
    node_size(node_count) = lay3_size;
end

%% Plot graph
G = graph(s, t);
fig = figure;
ax = axes('Position',[0.01 0.35 0.98 0.3]);
p = plot(ax,G, 'Layout', 'layered', 'NodeLabel',{});
p.Marker = 'o';
p.MarkerSize = node_size;
p.EdgeAlpha = 1.0;
G.Nodes.value = Node_value;
G.Edges.value = Edge_value;
G.Nodes.NodeColors = G.Nodes.value;
G.Edges.EdgeColors = G.Edges.value;
p.NodeColor = G.Nodes.NodeColors;
p.EdgeColor = G.Edges.EdgeColors;
p.LineWidth = 0.001;
set(gca,'box','off');
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
axis off;
            
end

