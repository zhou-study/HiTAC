function  [labels, np_list, pointE] = HiTAC(data, K)

data(:,all(data==0, 1))=[];
n=size(data,1);
data=data_norm(data);%Normalization
[knn_neigh, knn_neigh0, np_list, kf] = KNEG_pro(data);

% 构建图矩阵
pointE  = build_full_adj_matrix(knn_neigh(:, 2:kf));

np_list = localmax_degree(pointE);
np_list_new = remove_mutual_knn(np_list, knn_neigh0);
initC = data(np_list_new,:);
% 根据kmeans找索引
[labels, ~] = kmeans(data,size(initC,1), "Start",initC);

[~, ~, labels] = unique(labels, 'stable');





labelsC = 1:length(unique(labels));
[C, ~, ~, ~, ~] = ...
    cluster_flow_from_E(pointE, labels, 'IgnoreNonPositive', true, 'ZeroDiagonal', true);
d = balance_strength(C);
% 
% 
d = 1./d;
% % 找邻居
[cluster_neigh, ~] = top2_neighbors_from_d(d);
% 
[labelsC1, ~] = first( ...
    cluster_neigh, d, labelsC, K, pointE, ...
    'ShowPlot', false, ...
    'Data', data, ...
    'LabelsIn', labels, ...
    'Pause', 0.3, ...
    'DrawLink', true);
% % 第二次合并,标题的标签存在问题
[labelsC2, ~] = second( ...
    labelsC1, cluster_neigh, d, K, pointE,...
    'ShowPlot', false, 'Data', data, 'LabelsIn', labels, ...
    'Pause', 0.2, 'CircleSize', 12, 'CircleLineWidth', 2, 'DrawLink', true);
% % 第三次合并
[labelsC3, ~] = third( ...
    cluster_neigh, d, labelsC2, K, pointE,...
    'ShowPlot', false, ...
    'Data', data, ...
    'LabelsIn', labels, ...
    'Pause', 0.0, ...
    'DrawLink', true);






% 第四次合并
[labelsC_out, ~, ~] = last(labelsC3, d, K, pointE,...
    'ShowPlot', false, 'Data', data, 'LabelsIn', labels, 'Pause', 0.2, ...
    'DrawLink', true, 'CircleSize', 14, 'CircleLineWidth', 1.5);


labels= update_labels_from_labelsC(labels, labelsC_out);



end % function



function labels_new = update_labels_from_labelsC(labels, labelsC)
    labels  = labels(:);
    labelsC = labelsC(:);

    labels_new = labels;   % 先拷贝一份

    unique_super = unique(labelsC);
    unique_super(unique_super <= 0) = [];

    for s = unique_super'
        member_clusters = find(labelsC == s);   % 这个超簇包含的原簇 ID
        rep_id = min(member_clusters);         % 选一个代表ID（这里取最小）

        idx_points = ismember(labels, member_clusters);
        labels_new(idx_points) = rep_id;       % 点的标签改成代表ID
    end

    % labels <= 0（噪声）保持原值
    % 已经在开头拷贝时保留了，不用再处理
end