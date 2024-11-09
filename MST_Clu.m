function  [Label]=MST_Clu(Crep,dist,eta,M,cl)
%% -------------------------------------------------------------------------
% Aim:
% Construct MST on core representative points and final clustering
% -------------------------------------------------------------------------
% Input:
% Crep: core representative points
% dist: density-adaptive distance matrix
% eta: cluster size parameter
% M: the number of clusters
% -------------------------------------------------------------------------
% Output:
% Label: the clustering result
% -------------------------------------------------------------------------

[N,~]=size(Crep);
if N > M
    %% Construct MST on core representative points
    dist(dist == inf) = 0;
    G = sparse(dist);
    G = graph(G);
    [ST,pred] = minspantree(G,'Type','forest','Method','sparse'); % 'sparse' means that the Kruskal's (rather than Prim's) method is used
    ST = sparse(ST.Edges.EndNodes(:,1),ST.Edges.EndNodes(:,2),ST.Edges.Weight,N,N);
    
    if any(isnan(pred))
        % % Note: there is problem for the output 'pred' when UG is unconnected for Matlab 2013; this problem has been modified by Matlab 2016 and higher versions
        error('a higher version of Matlab is needed (e.g., Matlab 2016 and later version than that)');
    end
    
    idx_roots = find(pred == 0);  % transform MST to in-tree-based MSF
    pred(idx_roots) = idx_roots;
    omega = zeros(1,N);           % weight of core representative points
    for i = 1:N
        omega(i) = max(ST(i,pred(i)),ST(pred(i),i));  % ST is not a symmetric matrix, and thus max is used to get the non-zero element
    end
    
    %% Cut the tree and update cluster labels
    [rep,Crep] = EdgeCutting(pred,omega,cl,eta,M); 
    visited=zeros(N,1);  % updata representative points
    for i = 1:N
        if rep(i) ~= i
            parent = i;
            round = 1;
            visited(round) = i;
            while rep(parent) ~= parent
                parent = rep(parent);
                round = round+1;
                visited(round) = parent;
            end
            rep(visited(1:round)) = parent;
        end
    end
    Label = zeros(N,1); % Preallocate the space for the cluster label vector; Label(i): cluster label of point i;
    Label(Crep) = 1:length(Crep); % first, assign cluster labels to the core respresentative points;
    Label = Label(rep); % then, assign cluster labels to normal points based on the core respresentative points they have reached
    Label = Label(cl); % then, get the cluster label of all points
    
else
    warning(' in this case, cutting the tree is not needed')
    Label = cl;
end
end

