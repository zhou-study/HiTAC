function [C, cluster_ids, out_deg_cluster, in_deg_cluster, cluster_sizes] = cluster_flow_from_E(E, labels, varargin)
%CLUSTER_FLOW_FROM_E 将样本级完全邻接矩阵 E 汇总为簇级有向矩阵 C
%
% 输入
%   E        : N×N double，样本级完全邻接矩阵（0/1 或权重皆可），E(i,j) 表示 i→j
%   labels   : N×1 int，样本簇标签；相同标签属于同一簇。可包含 0 或负数表示噪声/忽略
%
% Name-Value（可选）
%   'IgnoreNonPositive' : logical，是否忽略 labels<=0 的样本（默认 true）
%   'ZeroDiagonal'      : logical，是否将簇内自环 C(i,i) 置 0（默认 false）
%
% 输出
%   C               : C×C double，簇级有向矩阵；C(i,j) 为簇 i → 簇 j 的边数/权重和
%   cluster_ids     : C×1 int，对应 C 的簇原始标签值（可能非连续）
%   out_deg_cluster : C×1 double，各簇出度（对所有簇的边数/权重和）
%   in_deg_cluster  : C×1 double，各簇入度
%   cluster_sizes   : C×1 int，各簇包含的样本数
%
% 说明
%   若 E 为二值矩阵，C 就是边数统计；若 E 为权重矩阵，C 为权重求和。
%   采用指示矩阵聚合：C = S' * E * S，其中 S(n,c)=1 表示样本 n 属于簇 c。

    p = inputParser;
    addParameter(p, 'IgnoreNonPositive', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'ZeroDiagonal', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    labels = labels(:);
    N = size(E,1);
    assert(size(E,2) == N, 'E 必须是 N×N 方阵');
    assert(numel(labels) == N, 'labels 长度必须与 E 维度一致');

    % 选择有效样本（可选忽略<=0标签）
    if opt.IgnoreNonPositive
        valid_mask = labels > 0;
    else
        valid_mask = true(N,1);
    end

    E_sub = E(valid_mask, valid_mask);
    lab_sub = labels(valid_mask);

    % 唯一簇ID（保持原始标签值顺序）
    cluster_ids = unique(lab_sub, 'stable');
    Cn = numel(cluster_ids);

    % 构建指示矩阵 S: (#valid N) × (Cn)
    S = zeros(sum(valid_mask), Cn);
    for c = 1:Cn
        S(:, c) = (lab_sub == cluster_ids(c));
    end

    % 簇级有向矩阵聚合：C = S' * E * S
    C = S' * E_sub * S;

    % 是否去掉簇内自环
    if opt.ZeroDiagonal
        C(1:Cn+1:end) = 0;
    end

    % 簇级入/出度与簇大小
    out_deg_cluster = sum(C, 2);      % 行和：簇出度
    in_deg_cluster  = sum(C, 1)';     % 列和：簇入度
    cluster_sizes   = sum(S, 1)';     % 每簇样本数
end