function [labelsC_out, merges] = second(labelsC, cluster_neigh, d, k, pointE,varargin)
% 第二次合并：基于"第一↔第二近邻"规则，按 d 的单链距离从小到大继续合并。
% 仅更新簇级标签 labelsC（行向量输出）；点级 labels 不更改。
%
% Inputs
%   labelsC        : 1xC 或 Cx1，第一次合并后的簇级标签（相同值 = 已同一超簇）
%   cluster_neigh  : Cx2，[first, second]；0 表示无该近邻
%   d              : CxC，原簇间距离矩阵
%   k              : 目标真实簇数
%
% Name-Value（可选，仅影响可视化，不改变逻辑）
%   'ShowPlot'       (false) : 是否画每步"合并前"图
%   'Data'           ([])    : n×2(+) 样本二维坐标
%   'LabelsIn'       ([])    : n×1 样本→原簇编号（1..C；若 0..C-1 会自动 +1）
%   'Pause'          (0.2)   : 每步暂停秒数
%   'CircleSize'     (10)    : 质心空心圆大小
%   'CircleLineWidth'(2)     : 质心圆线宽
%   'DrawLink'       (true)  : 是否画质心虚线
%
% Outputs
%   labelsC_out  : 1xC，更新后的簇级标签（与首次合并的"超簇ID体系"一致）
%   merges       : struct 数组（iter, from_super, to_super, dist）
%   finalK       : 最终簇数

    % ---------- 可视化参数 ----------
    p = inputParser;
    p.addParameter('ShowPlot', false, @(x)islogical(x) || isnumeric(x));
    p.addParameter('Data', [], @(x)isnumeric(x) && (isempty(x) || size(x,2)>=2));
    p.addParameter('LabelsIn', [], @(x)isnumeric(x) || islogical(x));
    p.addParameter('Pause', 0.2, @(x)isnumeric(x) && isscalar(x) && x>=0);
    p.addParameter('CircleSize', 10, @(x)isnumeric(x) && isscalar(x) && x>0);
    p.addParameter('CircleLineWidth', 2, @(x)isnumeric(x) && isscalar(x) && x>0);
    p.addParameter('DrawLink', true, @(x)islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opt = p.Results;

    % ---------- 基本检查 ----------
    labelsC = double(labelsC(:));
    C = numel(labelsC);
    if size(cluster_neigh,2) < 2, error('cluster_neigh 需至少两列 [first,second]'); end
    assert(all(size(d) == [C, C]), 'd 必须是 CxC');
    % 先设为恒等：每簇的父亲是自己
    parent = (1:C)';

    % 把"labelsC 值相等"的簇合并到同一个根（取该组最小索引为代表）
    pos_mask = isfinite(labelsC) & (labelsC > 0);    % 仅把 >0 的有效标签视为同簇先验
    uvals = unique(labelsC(pos_mask));
    for t = 1:numel(uvals)
        v = uvals(t);
        idx = find(labelsC == v);        % 该标签值对应的所有簇索引
        if numel(idx) <= 1, continue; end
        rep = min(idx);                  % 代表：最小索引（稳定且符合"并到较小ID"偏好）
        for ii = idx(:)'
            parent(ii) = rep;
        end
    end
    % 当前根数
    cur_roots = arrayfun(@(x) findroot(parent, x), (1:C)');
    curK = numel(unique(cur_roots));

    
    first  = double(cluster_neigh(:,1));
    second = double(cluster_neigh(:,2));
    targetK = max(1, min(k, numel(unique(labelsC))));

    % ---------- 1) 找"第一↔第二近邻"的原簇对 (i<j) ----------
    pairs_raw = [];
    for i = 1:C
        j = first(i);
        if j>0 && j<=C && second(j)==i && i<j
            pairs_raw = [pairs_raw; i, j]; %#ok<AGROW>
        end
        j = second(i);
        if j>0 && j<=C && first(j)==i && i<j
            pairs_raw = [pairs_raw; i, j]; %#ok<AGROW>
        end
    end
    if isempty(pairs_raw)
        labelsC_out = labelsC.';   % 行向量输出
        merges = struct('iter',{},'from_super',{},'to_super',{},'dist',{});
        finalK = numel(unique(labelsC));
        return;
    end
    % 去重并按距离升序排序
    pairs = unique(pairs_raw, 'rows', 'stable');
    cand_dist = arrayfun(@(r) d(pairs(r,1), pairs(r,2)), 1:size(pairs,1))';
    [cand_dist, ord ] = sort(cand_dist, 'ascend');
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

        % 画底图（尽量用你的 plot_clusters）
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

        % 点击显示样本序号
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
    % 规则：把较大的根并到较小的根（稳定代表 ID）
    if ra < rb
        parent(rb) = ra;
    else
        parent(ra) = rb;
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function txt = localUpdateFcn(evt, X)
    pos = evt.Position;
    idx = evt.DataIndex;  % 当前点在绘图数据中的索引
    % 查找该点在原始 Data 矩阵中的行号（若重复坐标则可能有多个）
    [~, data_idx] = ismember(pos, X, 'rows');
    if data_idx == 0
        data_idx = idx; % 如果没匹配到，就用绘图索引
    end
    txt = {sprintf('data行号 = %d', data_idx), ...
           sprintf('x = %.5g', pos(1)), ...
           sprintf('y = %.5g', pos(2))};
end

