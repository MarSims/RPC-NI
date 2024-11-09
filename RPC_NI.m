function [label] = RPC_NI(data, M, eta)
%% -------------------------------------------------------------------------
% Aim:
% The matlab code of the RPC-NI algorithm
% -------------------------------------------------------------------------
% Input:
% data: the dataset
% M: the number of clusters
% eta: cluster size parameter
% -------------------------------------------------------------------------
% Output:
% label: the clustering result
% -------------------------------------------------------------------------

%% Step 1: Determing core representative points
% disp('Determining core representative points');
[index,supk,Rho,rep,Crep,cl] = Rep_Point_Searching(data);
nCrep = length(Crep);
% disp('Number of core representation points：');disp(nCrep);

cdata=cell(1,nCrep);    % 保存每个簇中都有哪些点
cdataexp=cell(1,nCrep); % 保存每个簇中的点及每个点中的k近邻
nc=zeros(1,nCrep);      % 保存属于某个核心点的点数

for i=1:nCrep
    cdata{1,i}=find(rep==Crep(i))';
    [~,csize]=size(cdata{1,i});
    nc(i)=csize;
    temp=index(cdata{1,i},2:supk+1);
    if csize~=1
    cdataexp{1,i}=union(cdata{1,i},temp)';
    else
    temp2=union(temp,temp);
    cdataexp{1,i}=union(cdata{1,i},temp2);
    end
end

%% Step 2：Density-adaptive distance
% disp('Measuring density-adaptive distance matrix');
disCrep = pdist2(data(Crep,:),data(Crep,:));  % Euclidean distances between cores
CRho = Rho(Crep);
for i=1:nCrep
    for j=i+1:nCrep
        inset=intersect(cdataexp{1,i},cdataexp{1,j});    % common neighbors
        [~,numinset]=size(inset);
        if numinset==0
            disCrep(i,j) = exp(disCrep(i,j))-1;
        else
            sumrho = sum(Rho(inset));
            rhoCN = sumrho/numinset;
            DCC = (sqrt(CRho(i) * rhoCN) + sqrt(CRho(j) * rhoCN) + sqrt(CRho(i) * CRho(j))) / (CRho(i) + CRho(j) + rhoCN);
            disCrep(i,j) = exp(-DCC) * disCrep(i,j) / (sumrho);
        end
        disCrep(j,i) = disCrep(i,j);
    end
end

%% Step 3: MST-basd clustering
[label] = MST_Clu(data(Crep,:),disCrep,eta,M,cl);

end
