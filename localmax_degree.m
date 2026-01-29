function idx_peaks = localmax_degree(A)
% 从邻接矩阵 A 中找出 "度数局部极大" 的节点
% A(i,j)=1 表示 i→j
% 返回 idx_peaks：这些节点的索引

    % -------- 1. 基本统计 --------
    dout = full(sum(A, 2));      % 出度
    din  = full(sum(A, 1))';     % 入度
    deg  = dout + din;           % 总度（你定义的指标）
    N = size(A, 1);

    % -------- 2. 邻接定义（无向化）--------
    % 我们认为只要有一条边，就互为邻居
    Adj = (A > 0) | (A' > 0);
    Adj = Adj - diag(diag(Adj)); % 去掉自环

    % -------- 3. 判断局部极大 --------
    is_peak = false(N, 1);
    for i = 1:N
        nbr = find(Adj(i, :)); % 所有邻居
        if isempty(nbr)
            is_peak(i) = false; % 孤立点算峰
        else
            if all(deg(i) >= deg(nbr))
                is_peak(i) = true;
            end
        end
    end

    idx_peaks = find(is_peak);
end
