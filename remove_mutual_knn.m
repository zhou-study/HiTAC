function np_list_new = remove_mutual_knn(np_list, knn_neigh0)

keep = true(size(np_list));   % 标记是否保留

for a = 1:length(np_list)
    if ~keep(a)
        continue;
    end

    i = np_list(a);

    for b = a+1:length(np_list)
        if ~keep(b)
            continue;
        end

        j = np_list(b);

        % 判断是否互为近邻
        is_i_neighbor_j = any(knn_neigh0(i,:) == j);
        is_j_neighbor_i = any(knn_neigh0(j,:) == i);

        if is_i_neighbor_j && is_j_neighbor_i
            keep(b) = false;   % 删除后出现的那个
        end
    end
end

np_list_new = np_list(keep);

end
