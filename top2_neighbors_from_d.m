% 取近邻（排除 d==0）
function [cluster_neigh, cluster_dist] = top2_neighbors_from_d(d, tol)
%TOP2_NEIGHBORS_FROM_D  从距离矩阵 d 取每个簇的第一、第二近邻（排除0距离）
% 输入:
%   d   : CxC 距离矩阵（数值越小越近）。可非对称
%   tol : 可选，判断"等于0"的容差（默认 0，若有数值误差可设为1e-12等）
% 输出:
%   cluster_neigh : Cx2 int，邻居簇的索引（第一近邻、第二近邻）
%   cluster_dist  : Cx2 double，对应距离（若无则为 Inf）
%
% 规则:
%   - 忽略自身：对角线 = Inf
%   - 排除距离为0（或|d|<=tol）的条目：视为 Inf
%   - 若可用邻居不足2个，缺失位以 0（索引）与 Inf（距离）填充
%   - 并列(tie)时按索引小者优先（MATLAB排序稳定，保持该行为）

    if nargin < 2, tol = 0; end

    d = double(d);
    C = size(d,1);
    assert(size(d,2) == C, 'd 必须是 C×C 方阵');

    % 工作副本
    d_ = d;

    % 忽略自身
    d_(1:C+1:end) = Inf;

    % 排除等于0（或近似0）的距离
    if tol > 0
        mask_zero = abs(d_) <= tol;
    else
        mask_zero = (d_ == 0);
    end
    d_(mask_zero) = Inf;

    % 行内升序排序
    [vals, idx] = sort(d_, 2, 'ascend');

    % 取前两列（第一、第二近邻）
    if C >= 2
        take = min(2, max(0, C-1));   % 防御性处理
        cluster_neigh = idx(:, 1:take);
        cluster_dist  = vals(:, 1:take);
    else
        cluster_neigh = zeros(C, 2);
        cluster_dist  = Inf(C, 2);
        return;
    end

    % 若只有1个可用邻居，补齐到2列
    if size(cluster_neigh, 2) < 2
        cluster_neigh(:, end+1) = 0;
        cluster_dist(:,  end+1) = Inf;
    end

    % 没有近邻的位置（Inf）对应的索引置 0
    mask_no_neighbor = ~isfinite(cluster_dist);
    cluster_neigh(mask_no_neighbor) = 0;
end
