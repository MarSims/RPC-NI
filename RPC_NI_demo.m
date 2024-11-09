function RPC_NI_demo
% -------------------------------------------------------------------------
% Aim:
% The matlab code of "Representative Point-Based Clustering with Neighborhood
% Information for Complex Data Structures"
% -------------------------------------------------------------------------
% Input:
% data: the data set
% eta: cluster size parameter
% -------------------------------------------------------------------------
% Output:
% abel: the clustering result of local density peaks
% -------------------------------------------------------------------------
% Written by Zhongju Shang, Haowei Wang
% Novermber 2024

%% load dataset
 % clc; clear all;
data = textread('Jain.txt');
[N,dim]=size(data);
Tclass = data(:,dim);  % True labels
nClu = length(unique(Tclass));  % Number of clusters
data = data(:,1:dim-1); 

%% parameter
eta = 0.02*N;       % For most cases (eta controls the minimum cluster size); 
% eta = N/(2*nClu); % For equal size datasets, e.g., UCI datasets

%% clustering
[Label] = RPC_NI(data, nClu, eta); %% nC: number of clusters;
numClust = length(unique(Label));
disp('Number of Final Cluster');disp(numClust);

%% results evaluation
nmi = compute_nmi (Label,Tclass);
disp('NMI:');disp(nmi);

%% plot
figure;
clr = hsv(nClu*1);
gscatter(data(:,1),data(:,2),Tclass, clr, '.', 12);

end
