% 生成点级邻接矩阵E
function A = build_full_adj_matrix(knn_neigh)
%BUILD_FULL_ADJ_MATRIX 根据出度邻居矩阵 knn_neigh 生成有向邻接矩阵 A（非稀疏）
%
% 输入:
%   knn_neigh : N×K 整数矩阵
%                每行表示节点 i 的出度邻居（i → knn_neigh(i, j)）
%                若包含 0 或 i 自身，将被忽略。
%
% 输出:
%   A : N×N double 矩阵（完全矩阵形式）
%       A(i,j) = 1 表示 i → j 存在一条有向边，否则为 0。

    [N, K] = size(knn_neigh);
    A = zeros(N, N);   % 初始化为全零矩阵

    for i = 1:N
        neigh = knn_neigh(i, :);
        % 去除无效索引（0 或自身）
        neigh = neigh(neigh > 0 & neigh ~= i);
        A(i, neigh) = 1;
    end
end