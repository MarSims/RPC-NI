function  [cl]=MST_CLU(Crep,dist,nc,minsize,clu_num)
%% -------------------------------------------------------------------------
%Aim:
% Construct MST on Creps and cluster local densitypeaks
% -------------------------------------------------------------------------
% Input:
% A: core reps
% dist: density-adaptive distance matrix
% minsize: Minsize
% clu_num: the number of clusters
% -------------------------------------------------------------------------
% Output:
% cl: the clustering result of core representative points
% -------------------------------------------------------------------------

%% Construct MST on core reps
[N,~]=size(Crep);
dist(dist == inf) = 0;
G = sparse(dist);
G = graph(G);
[ST,pred] = minspantree(G,'Type','forest','Method','sparse'); % 'sparse' means that the Kruskal's method is used
ST = sparse(ST.Edges.EndNodes(:,1),ST.Edges.EndNodes(:,2),ST.Edges.Weight,N,N);

%========== Plot the minimum spanning tree ===========
% figure(1);
% % plot(Crep(:,1),Crep(:,2),'.');
% hold on;
% for i=1:N
%     j=pred(i);
%     if j~=0&&~isnan(j)
%         x=[Crep(i,1),Crep(j,1)];
%         y=[Crep(i,2),Crep(j,2)];
%         plot(x,y,'linewidth',1.5,'color','g');
%         hold on;
%     end
% end
% hold off;
%======================================================

% Determine optimal number of clusters
[xlabel,ylabel] = find(ST);
edges = full(ST(sub2ind(size(ST), xlabel, ylabel)));
[~,sortedind] = sort(edges,'descend');
ST2 = full(ST);

%=========== Plot the minimum spanning tree ============
%  figure(1);
%  plot(Crep(:,1),Crep(:,2),'.');
%  hold on;
% for i=1:N
%     j=pred(i);
%     if j~=0
%         x=[Crep(i,1),Crep(j,1)];
%         y=[Crep(i,2),Crep(j,2)];
%         plot(x,y);
%         % text(Crep(j,1),Crep(j,2),num2str(edges(j)));
%         hold on;
%     end
% end
%=======================================================

%% Cut the forest
i=0;
k=1;
while k < clu_num
    % 确定要去除的边
    s1=0;s2=0;
    p = xlabel(sortedind(i+1));
    q = ylabel(sortedind(i+1));
    while s1<minsize || s2<minsize
        s1=0;s2=0;
        visited = zeros(N,1);
        i = i+1;
        t = sortedind(i);
        xl(i) = xlabel(t);
        yl(i) = ylabel(t);
        p=xl(i);q=yl(i);

        queue=zeros(N,1);
        front=1;
        rear=1;
        queue(rear)=p;
        rear=rear+1;
        while front ~= rear
            temp=queue(front);
            s1=s1+nc(temp);
            visited(temp)=1;
            front=front+1;
            for j=1:N
                if (ST2(temp,j)~=0||ST2(j,temp)~=0)&&j~=q&&visited(j)==0&&iscontain2(queue,front,rear,j)==0
                    queue(rear)=j;
                    rear=rear+1;
                end
            end
        end
        queue=zeros(N,1);
        front=1;
        rear=1;
        queue(rear)=q;
        rear=rear+1;
        while front~=rear
            temp=queue(front);
            s2=s2+nc(temp);
            visited(temp)=1;
            front=front+1;
            for j=1:N
                if (ST2(temp,j)~=0||ST2(j,temp)~=0)&&j~=p&&visited(j)==0&&iscontain2(queue,front,rear,j)==0
                    queue(rear)=j;
                    rear=rear+1;
                end
            end
        end  
    end

    ST2(p,q)=0;
    ST2(q,p)=0;
    k=k+1; 
end

cl=zeros(N,1);
   ncl=0;
   sumedge=zeros(k,1);
   for i=1:N
       if cl(i)==0
           ncl=ncl+1;
           queue=zeros(N,1);
           front=1;
           rear=1;
           queue(rear)=i;
           rear=rear+1;
           sumedge(ncl)=0;
           while front~=rear
               p=queue(front);
               front=front+1;
               cl(p)=ncl;
               for j=1:N
                   if (ST2(p,j)~=0||ST2(j,p)~=0)&&cl(j)==0&&iscontain(queue,j)==0
                       queue(rear)=j;
                       rear=rear+1;
                   end
               end
           end
       end
   end
end




