function [pr2,rs2] = Edge_Cutting(pred,weightE,c,eta,M)
%============================================================
% Aim:
% MST-based clustering
% ----------------------------------------------------------
% Input:
% pred: the parent (core) point of (core) point i
% weightE: the edge weight vector
% c: initial cluster label vector
% eta: cluster size parameter
% M: the number of clusters
% ----------------------------------------------------------
% Output:
% pr: the updated representative vector after edge cutting.
% rs: the core representative vector. 
%=============================================================

%% determine the initial weight of each supernode (see literature [23] for details)
len = length(weightE); N = length(c);
NW =zeros(1,len); % weightCr(i) will store the number of points in cluster i.
for i = 1:N
    NW(c(i)) = NW(c(i)) + 1;
end

inNei = cell(1,len); % note: each cell is not initialized by the array size; this is a problem for speed as the elements are gradually added; but this is not a problem for c++ using advanced data structure; For Matlab, one can also consider to use a vector to count the number of elements of each cell first, based on which each cell can be initialized; 
for i = 1:len
    if pred(i) ~= i
        inNei{1,pred(i)}=[inNei{1,pred(i)},i];
    end
end

i = 1;
if size(pred,1)>size(pred,2)
    pred = pred';
end
roots = find(pred == (1:len));

Ind = zeros(len,1);
ct = 0;
for j = 1:length(roots)       
    ct = ct + 1;
    Ind(ct) = roots(j); 
    while i<=ct
        temp = inNei{1,Ind(i)};
        if ~isempty(temp)     
            for m = 1:length(temp)                
                ct = ct + 1;
                Ind(ct) = temp(m);
            end
        end
        i = i + 1;
    end
end
                
% update the weight of each core resentative points on the in-tree-based forest
for j = len:-1:1  
    i = Ind(j);
    if pred(i)~=i
        NW(pred(i))=NW(pred(i))+NW(i);
    end
end

%% cut the edges according to cluster number and MinSize
[~,ind]=sort(weightE,'descend');
cutE = 0;       % number of cut edge
m = 1;          % number of clusters in MST
num_roots_initial = length(roots);
num_of_edges_in_graph = len - num_roots_initial;
if M >= num_roots_initial
    num_of_edges_required_to_remove = M - num_roots_initial; 
else
    num_of_edges_required_to_remove = 0;
    warning('there could exist over-partitioning problem; it is suggest to increase the value of parameter k or increase cluster number');
end
passed_node = zeros(len,1);

% disp('check the edges one by one in decreasing order of edge weight...')
while cutE ~= num_of_edges_required_to_remove && m <= num_of_edges_in_graph
    start_node = ind(m);
    end_node = pred(start_node);
    if NW(start_node) > eta
        ct = 1; 
        passed_node(ct) = end_node;  
        while end_node ~= pred(end_node)
            end_node = pred(end_node);  
            ct = ct + 1;
            passed_node(ct) = end_node;
        end
        root_node_reached = end_node;
        if NW(root_node_reached)-NW(start_node) > eta
            NW(passed_node(1:ct)) = NW(passed_node(1:ct)) - NW(start_node);
            pred(start_node) = start_node;
            cutE = cutE + 1; 
        end
    end
    m = m + 1;
end

if cutE < num_of_edges_required_to_remove
    warning('The value of eta is too large, preventing further edge cuts. Consider reducing eta.');
end

pr2 = pred;
rs2 = find(pr2 == (1:len)); 

end
