function [labelsC_out, merges, finalK] = last(labelsC, d, k, pointE,varargin)
%MERGE_BY_D_UNTIL_K
% 根据距离矩阵 d，从最小的簇对开始，依次在 labelsC 中标记为同一簇（不改点级 labels），
% 已经在同一超簇则跳过；当不同标签数 = k 时停止。
%
% Inputs
%   labelsC : 1×C 或 C×1，簇级"同簇关系"标记（同值表示同一超簇；初始通常为 1..C）
%   d       : C×C，簇间距离矩阵（越小越相似，d(i,i) 可为 0；允许 Inf/NaN 将被跳过）
%   k       : 目标超簇数（不同标签个数）
%
% Name-Value 可选参数（与可视化相关）
%   'ShowPlot'       : 是否逐步可视化（默认 false）
%   'Data'           : n×2(+) 样本坐标（仅用于绘图）
%   'LabelsIn'       : n×1 样本→原簇编号（1..C；若 0..C-1 则自动 +1）
%   'Pause'          : 每步暂停秒数（默认 0.2）
%   'DrawLink'       : 是否在两簇质心间画虚线（默认 true）
%   'CircleSize'     : 质心圈大小（默认 14）
%   'CircleLineWidth': 质心圈线宽（默认 1.5）
%
% Outputs
%   labelsC_out : 1×C，最终的簇级"同簇关系"标记（每个簇被映射为其超簇代表ID = 超簇内最小索引）
%   merges      : 结构体数组，记录每步（候选对、距离、合并前后簇数等）
%   finalK      : 停止时的不同标签个数
%
% 说明
% - "不因为两个簇标记为同一个簇就更新 labels（点级）"，因此本函数只改变 labelsC，不改变点级 labels。
% - 如果 1 号簇和 2 号簇合并，则最终 labelsC(2) 会等于 1（因为我们把每个超簇的代表定为超簇内最小索引）。
%

%% -------- 参数检查与读取 --------
if isrow(labelsC), labelsC = labelsC(:); end
C = numel(labelsC);
if ~ismatrix(d) || size(d,1)~=C || size(d,2)~=C
    error('d 大小应为 C×C，且与 labelsC 一致。');
end
p = inputParser;
p.addParameter('ShowPlot', false, @(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('Data', [], @(x)isnumeric(x)&&~isempty(x));
p.addParameter('LabelsIn', [], @(x)isnumeric(x)||isempty(x));
p.addParameter('Pause', 0.2, @isscalar);
p.addParameter('DrawLink', true, @(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('CircleSize', 14, @isscalar);
p.addParameter('CircleLineWidth', 1.5, @isscalar);
p.parse(varargin{:});
opt = p.Results;

%% -------- 用并查集承载"同簇关系" --------
% 初始 parent：每个簇一个根
parent = (1:C).';
rankUF = zeros(C,1,'uint8');   % 可选：按秩合并
minrep = (1:C).';              % 维护每个根代表的"最小索引"，确保最终 labelsC 映射为超簇内最小ID

% 若 labelsC 中已有"同值=同簇"，先把它们并起来
[~,~,gid] = unique(labelsC(:),'stable');
for g = unique(gid).'
    members = find(gid==g);
    if numel(members)>=2
        r = members(1);
        for j = 2:numel(members)
            [parent, rankUF, minrep] = uf_union(parent, rankUF, minrep, r, members(j));
            r = uf_find(parent, r); % 更新当前根
        end
    end
end

% 当前不同超簇数
curK = numel(unique(arrayfun(@(c) uf_find(parent,c), 1:C)));
if curK <= k
    labelsC_out = compress_labelsC(parent, minrep, C);
    finalK = curK;
    merges = struct('iter',{},'a',{},'b',{},'dist',{},'beforeK',{},'afterK',{});
    return;
end

%% -------- 构造候选簇对（按 d 从小到大）--------
D = d;

% 忽略下三角和对角线
D(tril(true(C))) = Inf;

% 忽略 NaN
D(isnan(D)) = Inf;

% 忽略距离为 0（或近似 0）的情况
tol = 1e-12;                % 容差阈值，可按需要调整
D(abs(D) <= tol) = Inf;

% 提取候选对：仅上三角、有限值
mask  = triu(true(C), 1) & isfinite(D);
lin   = find(mask);          % 线性索引
vals  = D(lin);              % 对应距离

% 按距离从小到大排序
[sorted_dist, order] = sort(vals, 'ascend');
[i, j] = ind2sub([C, C], lin(order));

% 输出候选簇对及距离
pairs = [i, j];
cand_dist = sorted_dist;

%% -------- 逐对尝试合并（只改 labelsC/并查集；不改点级 labels）--------
merges = struct('iter',{},'a',{},'b',{},'dist',{},'beforeK',{},'afterK',{});
t = 0;
for i = 1:size(pairs,1)
    a_id = pairs(i,1);
    b_id = pairs(i,2);
    ra = uf_find(parent, a_id);
    rb = uf_find(parent, b_id);

    % 已在同一超簇 -> 跳过
    if ra == rb, continue; end

    % ---------- 合并前可视化（不改点级 labels） ----------
    if opt.ShowPlot && ~isempty(opt.Data) && ~isempty(opt.LabelsIn)
    
        X = opt.Data(:,1:2);
        LabIn = double(opt.LabelsIn(:));

        % 若点级标签为 0..C-1，则自动加 1 到 1..C
        if ~isempty(LabIn)
            if min(LabIn)==0 && max(LabIn)==C-1, LabIn = LabIn + 1; end
        end
        valid_pts = (LabIn >= 1) & (LabIn <= C);

        % 每个原簇(1..C)当前的根超簇ID
        root_of_super = zeros(C,1);
        for c0 = 1:C
            root_of_super(c0) = uf_find(parent, c0);
        end
        labels_plot_disp = zeros(size(LabIn));
        labels_plot_disp(valid_pts) = root_of_super(LabIn(valid_pts));

        % % 画底图（优先用你的 plot_clusters）
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

        t_disp = sum([merges.iter]==1) + 1; % 仅用于标题显示的递增序号
        title(sprintf('Iter %d (pre-merge): C%d <-> C%d | dist=%.6g | K\\rightarrow%d', ...
              t_disp, a_id, b_id, cand_dist(i), curK-1), 'Interpreter','tex');

        % 点击显示样本信息
        dcm = datacursormode(gcf);
        set(dcm,'Enable','on','UpdateFcn',@(obj,evt) localUpdateFcn(evt, X));

        drawnow;
        if opt.Pause>0, pause(opt.Pause); end
    end
    % ---------- 可视化结束 ----------

    % 真正执行合并（只更新并查集）
    beforeK = curK;
    [parent, rankUF, minrep] = uf_union(parent, rankUF, minrep, ra, rb);

    % 更新当前不同超簇数
    curK = numel(unique(arrayfun(@(c) uf_find(parent,c), 1:C)));

    % 记录日志
    t = t + 1;
    merges(t).iter    = 1;          % 每条记录代表一次有效合并
    merges(t).a       = a_id;
    merges(t).b       = b_id;
    merges(t).dist    = cand_dist(i);
    merges(t).beforeK = beforeK;
    merges(t).afterK  = curK;

    % 达到目标 k -> 停止
    if curK == k
        break;
    end
end

%% -------- 输出 labelsC（每簇映射到其超簇代表：超簇内最小索引）--------
labelsC_out = compress_labelsC(parent, minrep, C);
finalK = curK;

end

%% ======================= 辅助函数 =======================

function r = uf_find(parent, x)
% 带路径压缩
while parent(x) ~= x
    parent(x) = parent(parent(x));
    x = parent(x);
end
r = x;
end

function [parent, rankUF, minrep] = uf_union(parent, rankUF, minrep, x, y)
% 按秩合并，并维护每个根的"超簇代表最小索引"
rx = uf_find(parent, x);
ry = uf_find(parent, y);
if rx == ry, return; end
% 统一让更高秩作为根
if rankUF(rx) < rankUF(ry)
    parent(rx) = ry;
    minrep(ry) = min(minrep(ry), minrep(rx));
elseif rankUF(rx) > rankUF(ry)
    parent(ry) = rx;
    minrep(rx) = min(minrep(rx), minrep(ry));
else
    parent(ry) = rx;
    rankUF(rx) = rankUF(rx) + 1;
    minrep(rx) = min(minrep(rx), minrep(ry));
end
end

function labelsC_out = compress_labelsC(parent, minrep, C)
% 将每个簇 i 映射到其根 r 的代表（超簇内最小索引）
labelsC_out = zeros(C,1);
for i = 1:C
    r = uf_find(parent, i);
    labelsC_out(i) = minrep(r); % 满足"若 1 与 2 合并，则 labelsC(2)=1"
end
labelsC_out = labelsC_out(:).';
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
