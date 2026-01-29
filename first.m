function [labelsC_out, info] = first(cluster_neigh, d, labelsC, k, pointE,varargin)
% MERGE_BY_MUTUAL_FIRST_UNTIL_K
% 仅依据"互为第1近邻"的簇对，按距离矩阵 d 从小到大依次合并（只更新簇级 labelsC，不改动点级 labels，也不更新邻接）
%
% Inputs
%   cluster_neigh : C×K int，簇的邻接表（升序；第1列为"最近邻(第1近邻)"，1-based）
%   d             : C×C double，簇间距离（对称；d(i,i)=0）
%   labelsC       : C×1 int，"同值同簇"的簇标签（不是父指针；相同值表示已属同一超簇）
%   k             : 目标簇数（当根簇数<=k 或无候选对时停止）
%
% Name-Value（可选）
%   'ShowPlot'        : logical，逐步可视化（默认 false）
%   'Data'            : n×2/≥2 double，原始数据坐标（可视化用）
%   'LabelsIn'        : n×1 int，点级簇标签（仅映射显示，不会被修改）
%   'Pause'           : double，合并前暂停秒数（默认 0）
%   'DrawLink'        : logical，是否画两簇质心连线（默认 true）
%   'CircleSize'      : double，质心圆大小（默认 10）
%   'CircleLineWidth' : double，质心圆线宽（默认 1.5）
%
% Outputs
%   labelsC_out : C×1 int，合并后的簇级标签（用"根簇索引"表示；同簇同值）
%   info        : struct，过程信息（候选对、距离、合并日志、停止原因等）

% -------- 解析选项 --------
ip = inputParser;
ip.addParameter('ShowPlot', false);
ip.addParameter('Data', []);
ip.addParameter('LabelsIn', []);
ip.addParameter('Pause', 0);
ip.addParameter('DrawLink', true);
ip.addParameter('CircleSize', 10);
ip.addParameter('CircleLineWidth', 1.5);
ip.parse(varargin{:});
opt = ip.Results;

% -------- 基础检查 --------
labelsC = double(labelsC(:));
C = numel(labelsC);

if size(cluster_neigh,2) < 2, error('cluster_neigh 需至少两列 [first,second]'); end
assert(all(size(d) == [C, C]), 'd 必须是 CxC');

% =========================================================
% 1) 初始化并查集（注意：labelsC 是"同值同簇"，不是父指针）
% =========================================================
parent = (1:C)';                             % 先设为恒等
pos_mask = isfinite(labelsC) & (labelsC > 0);% 仅把 >0 的标签视为同簇先验
uvals = unique(labelsC(pos_mask));
for t = 1:numel(uvals)
    v = uvals(t);
    idx = find(labelsC == v);                % 该标签值对应的所有簇索引
    if numel(idx) <= 1, continue; end
    rep = min(idx);                          % 代表：最小索引（稳定）
    for ii = idx(:)'
        parent(ii) = rep;
    end
end

% 当前根数
cur_roots = arrayfun(@(x) findroot(parent, x), (1:C)');
curK = numel(unique(cur_roots));

% =========================================================
% 2) 构造"互为第1近邻"的候选合并对
% =========================================================
if size(cluster_neigh,2) < 1
    labelsC_out = arrayfun(@(i) findroot(parent, i), (1:C)');
    info = struct('cand_pairs_super', [], 'cand_dist', [], ...
              'merges', [], 'stoppedReason', 'no_first_neighbor_column', ...
              'finalK', curK);
    return;
end

firstN = cluster_neigh(:, 1);                        % 第1近邻
valid  = firstN>=1 & firstN<=C & firstN ~= (1:C)';   % 合法且不是自己
pairs = [];
for i = 1:C
    if ~valid(i), continue; end
    j = firstN(i);
    if j<1 || j>C || i==j, continue; end
    % 互为第1近邻
    if cluster_neigh(j,1) == i
        a = min(i,j); b = max(i,j);
        pairs = [pairs; a b]; %#ok<AGROW>
    end
end

if isempty(pairs)
    labelsC_out = arrayfun(@(i) findroot(parent, i), (1:C)');
    info = struct('cand_pairs_super', [], 'cand_dist', [], 'merges', [], ...
                  'stoppedReason', 'no_mutual_first_neighbors', ...
                  'finalK', curK);
    return;
end

% 去重并按距离升序排序
pairs = unique(pairs, 'rows', 'stable');
cand_dist = arrayfun(@(r) d(pairs(r,1), pairs(r,2)), 1:size(pairs,1))';
[ cand_dist, ord ] = sort(cand_dist, 'ascend');
cand_pairs_super = pairs(ord, :);

% =========================================================
% 3) 主循环：按距离从小到大尝试合并（仅更新簇级 parent）
% =========================================================
merges = struct('step', {}, 'from', {}, 'to', {}, 'dist', {}, 'curK_after', {});
targetK = k;

for t = 1:numel(cand_dist)
    if curK <= targetK, break; end

    a_id = cand_pairs_super(t,1);
    b_id = cand_pairs_super(t,2);

    ra = findroot(parent, a_id);
    rb = findroot(parent, b_id);
    if ra == rb
        continue; % 已经同簇
    end

    % ---------- 合并前可视化（不改点级 labels） ----------
    if opt.ShowPlot && ~isempty(opt.Data) && ~isempty(opt.LabelsIn)
   
        X = opt.Data(:,1:2);
        LabIn = double(opt.LabelsIn(:));

        % 若点级标签为 0..C-1，则自动加 1 到 1..C
        if ~isempty(LabIn)
            if min(LabIn)==0 && max(LabIn)==C-1, LabIn = LabIn + 1; end
        end
        valid_pts = (LabIn >= 1) & (LabIn <= C);

        % 每个"原簇ID（1..C）"当前的根超簇 ID（用于把点映射到当前超簇）
        root_of_super = zeros(C,1);
        for c0 = 1:C
            root_of_super(c0) = findroot(parent, c0);
        end
        labels_plot_disp = zeros(size(LabIn));
        labels_plot_disp(valid_pts) = root_of_super(LabIn(valid_pts));

        % 画底图（优先用你的 plot_clusters）
        % newFig = figure; %#ok<NASGU>
        try
            plot_clusters(X, labels_plot_disp, pointE); hold on;
            
        catch
            % 兜底绘制
            hold on; box on; axis equal;
            uLabs = unique(labels_plot_disp); uLabs(uLabs==0) = [];
            cmap = lines(max(1,numel(uLabs)));
            for ii = 1:numel(uLabs)
                mk = labels_plot_disp == uLabs(ii);
                scatter(X(mk,1), X(mk,2), 12, cmap(ii,:), 'filled', ...
                        'HitTest','on','PickableParts','all');
            end
        end

        % 把 labels_plot_disp 放入 Figure.UserData，便于 Data Cursor 读取
        set(gcf, 'UserData', struct('labels_plot_disp', labels_plot_disp));

        % 两个超簇（ra, rb）的样本索引与质心
        idx_a = find(opt.LabelsIn == a_id);
        idx_b = find(opt.LabelsIn == b_id);
        have_a = any(idx_a); have_b = any(idx_b);
        % if have_a
        %     ca = mean(X(idx_a,:),1);
        %     plot(ca(1), ca(2), 'ko', 'MarkerSize', opt.CircleSize, ...
        %          'LineWidth', opt.CircleLineWidth, 'MarkerFaceColor','none');
        % end
        % if have_b
        %     cb = mean(X(idx_b,:),1);
        %     plot(cb(1), cb(2), 'ko', 'MarkerSize', opt.CircleSize, ...
        %          'LineWidth', opt.CircleLineWidth, 'MarkerFaceColor','none');
        % end
        % if opt.DrawLink && have_a && have_b
        %     plot([ca(1) cb(1)], [ca(2) cb(2)], 'k--', 'LineWidth', 1.2);
        % end

        title(sprintf('Iter %d (pre-merge): C%d <-> C%d | dist=%.6g | K\\rightarrow%d', ...
              t, a_id, b_id, cand_dist(t), curK-1), 'Interpreter','tex');

        % 点击显示样本序号（以及 data 行号与当前超簇ID）
        dcm = datacursormode(gcf);
        set(dcm,'Enable','on','UpdateFcn',@(obj,evt) localUpdateFcn(evt, X));

        drawnow;
        if opt.Pause>0, pause(opt.Pause); end
    end
    % ---------- 可视化结束 ----------

    % ---- 执行合并（只改簇级 parent/根关系）----
    parent = union_roots(parent, ra, rb);

    % 更新当前 K（根数）
    cur_roots = arrayfun(@(x) findroot(parent, x), (1:C)');
    curK = numel(unique(cur_roots));

    % 记录日志
    merges(end+1).step = t; %#ok<AGROW>
    merges(end).from   = rb;
    merges(end).to     = ra;
    merges(end).dist   = cand_dist(t);
    merges(end).curK_after = curK;
end

% =========================================================
% 4) 输出与信息
% =========================================================
labelsC_out = arrayfun(@(i) findroot(parent, i), (1:C)'); % 用"根簇索引"表达，同簇同值
stoppedReason = ternary(curK <= targetK, 'reach_target_k', 'exhaust_pairs');

info = struct();
info.cand_pairs_super = cand_pairs_super;
info.cand_dist        = cand_dist;
info.merges           = merges;
info.finalK           = curK;
info.stoppedReason    = stoppedReason;

end % ====== 主函数结束 ======


% ====== 工具函数们 ======
function r = findroot(parent, i)
    r = i;
    while parent(r) ~= r
        parent(r) = parent(parent(r)); % 路径压缩
        r = parent(r);
    end
end

function parent = union_roots(parent, a, b)
    ra = findroot(parent, a);
    rb = findroot(parent, b);
    if ra == rb, return; end
    % 规则：把较大的根并到较小的根（代表ID稳定）
    if ra < rb
        parent(rb) = ra;
    else
        parent(ra) = rb;
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

% Data Cursor 回调：显示 data 行号 + 坐标 + 当前超簇ID（若可得）
function txt = localUpdateFcn(evt, X)
    pos = evt.Position;
    idx = evt.DataIndex;                 % 绘图数据中的索引
    [~, data_idx] = ismember(pos, X, 'rows');
    if data_idx == 0, data_idx = idx; end

    fig = ancestor(evt.Target, 'figure');
    clu = '';
    if isfield(fig.UserData, 'labels_plot_disp')
        lp = fig.UserData.labels_plot_disp;
        if data_idx >=1 && data_idx <= numel(lp)
            clu = sprintf('簇 = %d', lp(data_idx));
        end
    end
    if isempty(clu)
        txt = {sprintf('data行号 = %d', data_idx), ...
               sprintf('x = %.5g', pos(1)), ...
               sprintf('y = %.5g', pos(2))};
    else
        txt = {sprintf('data行号 = %d', data_idx), ...
               sprintf('x = %.5g', pos(1)), ...
               sprintf('y = %.5g', pos(2)), ...
               clu};
    end
end
