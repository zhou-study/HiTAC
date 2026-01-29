function [knn_neigh, knn_neigh0, np_list, kf] = KNEG_pro(data)

data(:,all(data==0, 1))=[];
n=size(data,1);
data=data_norm(data);%Normalization
k = floor(sqrt(n));
%% Construct KNEG
t1=clock;
[knn_neigh,knn_dist] = knnsearch(data,data,'k',k);
knn_neigh0=knn_neigh;%Backup result
knn_oppsite=ones(n,k);
for i=2:k
    for j=1:n
        % 找j的第i个邻居
        neigh=knn_neigh(j,i);
        if neigh~=0
            % 如果我有这个邻居，看看它是我的第几邻居
            loc=find(knn_neigh(neigh,:)==j);
            % 如果它不是我的邻居，或者loc大于i
            if isempty(loc) || loc>i
                % neigh的第loc个邻居为空
                knn_neigh(neigh,loc)=0;
                knn_dist(neigh,loc)=inf;
                knn_oppsite(j,i)=0;
            end
        end
    end
end

%% obtain peak density point
n_np_list=[];
is_np=[];
is_np2=[];
is_np_list=[];
np_list=[];
% last_flag=0;
for kk=1:k
    is_np=zeros(n,1);
    knn_oppsite2=~knn_oppsite;
    judge=sum(knn_neigh(:,1:kk).*knn_oppsite2(:,1:kk),2);
    is_np(judge==0)=1;
    is_np_list=[is_np_list,is_np];
    n_np=length(find(is_np==1));
    n_np_list=[n_np_list,n_np];
    if kk~=1
        if n_np_list(end-1)-n_np_list(end)<=n^0.5 %suspensive condition
            kf=kk;
            np_list=(find(is_np==1));
            break;
        end
    end
end
% np_list = find(is_np_list(:,2) == 1);
end