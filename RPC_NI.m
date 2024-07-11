function [cl] = RPC_NI(A, nClu, MinSize)
%% -------------------------------------------------------------------------
%Aim:
%The matlab code of the RPC-NI algorithm
% -------------------------------------------------------------------------
%Input:
%A: the data set
%clu_num:the number of clusters
% -------------------------------------------------------------------------
%Output:
%cl: the clustering result
%time: the running time of the RPC-NI
% -------------------------------------------------------------------------

%% Step1: Determing core representative points
[N,~]=size(A);
% disp('Determining local density peaks');
[index,supk,Rho,rep,Crep] = Rep_Point_Searching(A);

%% Step 2: Density-adaptive distance
% disp('Searching representative domain')
nCrep = length(Crep);
cdata=cell(1,nCrep); 
cdataexp=cell(1,nCrep); 
nc=zeros(1,nCrep); 

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

% disp('Computing density-adaptive distance');
disCrep = pdist2(A(Crep,:),A(Crep,:)); 
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


%% Step 3: MST-based clustering
% disp('Construct MST on local density peaks');
[cores_cl]=MST_CLU(A(Crep,:),disCrep,nc,MinSize,nClu);

cl=zeros(N,1);
for i=1:nCrep
    cl(Crep(i))=cores_cl(i);
end
for i=1:N
    cl(i)=cl(rep(i));
end
disp('Complete!');

end
