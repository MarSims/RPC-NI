function RPC_NI_demo
    %% load dataset
     % clc; clear all;
    data = textread('Jain.txt');
    [N,dim]=size(data);
    Tclass = data(:,dim);  % 真实类别
    nClu = length(unique(Tclass));  %类别数
    data = data(:,1:dim-1);   % 去除标签 
    
    %% parameter
    % MinSize = 0.02*N; % Note: parameter MinSize (i.e.,the minimal cluster size) is dependent on ratio; 
    MinSize = N/(2*nClu); % For equal size datasets, e.g., UCI datasets
    
    %% clustering
    [Label] = RPC_NI(data, nClu, MinSize); %% nC: number of clusters;
    
    %% result
    numClust = length(unique(Label));
    disp('Number of Final Cluster');disp(numClust);
    
    nmi = compute_nmi (Label,Tclass);
    disp('NMI:');disp(nmi);
    
    %% plot
    figure;
    clr = hsv(nClu*1);
    gscatter(data(:,1),data(:,2),Tclass, clr, '.', 12);

end