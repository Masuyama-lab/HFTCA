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
function Model2 = Generate_Model(Model)

% GENERATE_MODEL2  Change the order of model nodes in hue order.

Model2 = Model;

% lay1
[A1,B1] = size(Model2.Means);
for i = 1:B1
    for j = i+1:B1
        C = Model2.Means(:,i);
        D = Model2.Means(:,j);
        C = rgb2hsv(C.');
        D = rgb2hsv(D.');
        if C(1,1) > D(1,1)
            % Change Means & Child
            tmp_Means = Model2.Means(:,i);
            Model2.Means(:,i) = Model2.Means(:,j);
            Model2.Means(:,j) = tmp_Means;

            tmp_Child = Model2.Child(1,i);
            Model2.Child(1,i) = Model2.Child(1,j);
            Model2.Child(1,j) = tmp_Child;
        end
    end
end

% lay2
for child_idx = 1:B1
    if ~isempty(Model2.Child{1,child_idx})
        [A2,B2] = size(Model2.Child{1,child_idx}.Means);
        for i = 1:B2
            for j = i+1:B2
                C = Model2.Child{1,child_idx}.Means(:,i);
                D = Model2.Child{1,child_idx}.Means(:,j);
                C = rgb2hsv(C.');
                D = rgb2hsv(D.');
                if C(1,1) > D(1,1)
                    % Change Means & Child
                    tmp_Means = Model2.Child{1,child_idx}.Means(1,i);
                    Model2.Child{1,child_idx}.Means(1,i) = Model2.Child{1,child_idx}.Means(1,j);
                    Model2.Child{1,child_idx}.Means(1,j) = tmp_Means;

                    tmp_Child = Model2.Child{1,child_idx}.Child(1,i);
                    Model2.Child{1,child_idx}.Child(1,i) = Model2.Child{1,child_idx}.Child(1,j);
                    Model2.Child{1,child_idx}.Child(1,j) = tmp_Child;
                end
            end
        end
    end
end

% lay3
for child_idx = 1:B1
    if ~isempty(Model2.Child{1,child_idx})
        [A2,B2] = size(Model2.Child{1,child_idx}.Means);
        for child_child_idx = 1:B2
            if ~isempty(Model2.Child{1,child_idx}.Child{1,child_child_idx})
                [A3,B3] = size(Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means);
                for i = 1:B3
                    for j = i+1:B3
                        C = Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means(:,i);
                        D = Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means(:,j);
                        C = rgb2hsv(C.');
                        D = rgb2hsv(D.');
                        if C(1,1) > D(1,1)
                            % Change Means & Child
                            tmp_Means = Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means(1,i);
                            Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means(1,i) = Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means(1,j);
                            Model2.Child{1,child_idx}.Child{1,child_child_idx}.Means(1,j) = tmp_Means;

                            tmp_Child = Model2.Child{1,child_idx}.Child{1,child_child_idx}.Child(1,i);
                            Model2.Child{1,child_idx}.Child{1,child_child_idx}.Child(1,i) = Model2.Child{1,child_idx}.Child{1,child_child_idx}.Child(1,j);
                            Model2.Child{1,child_idx}.Child{1,child_child_idx}.Child(1,j) = tmp_Child;
                        end
                    end
                end
            end
        end
    end
end

end

