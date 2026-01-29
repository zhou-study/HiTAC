function plot_clusters(data, labels, pointE, base_cmap)
% data: n×2 坐标
% labels: n×1 整数簇编号（>0 为簇；0/负数表示未分配/噪声）
% pointE: n×n 邻接矩阵，非零表示 i -> j 存在边（可为0/1或带权重）
% base_cmap: 可选，string 数组/char 数组形式的十六进制颜色表

    figure; hold on; box on;

    n = size(data, 1); %#ok<NASGU>

    % -------- 0. 默认颜色表（如果没有传入） --------
    if nargin < 4 || isempty(base_cmap)
        base_cmap  = [
            "#FF0000";  % 1 红色
            "#FFA500";  % 2 橙色
            "#FFFF00";  % 3 黄色
           
    
            "#CCFF00";  % 4 柠檬绿
            "#00FF00";  % 5 绿色
            "#00FFFF";  % 6 青色
            "#007FFF";  % 7 蔚蓝
            "#6B8E23";  % 18 橄榄绿
            
            "#7F00FF";  % 14 紫罗兰 
         
            "#FF00FF";  % 10 洋红
            "#FFC0CB";  % 11 粉红
            "#DDA0DD";  % 12 李子紫
            "#800080";  % 13 紫色
            "#000080";  % 9 海军蓝
            "#00BFFF";  % 15 深天蓝
            "#40E0D0";  % 17 绿宝石
            "#7FFFD4";  % 8 蓝绿
            "#808000";  % 19 橄榄黄
            "#F5DEB3";  % 20 小麦色
            "#008000";  % 21 深绿
            "#8B4513";  % 22 马鞍棕
            "#E3F9FD";  % 23 莹白
            "#AA4C8F";  % 24 梅紫
            "#DF7163";  % 25 浅绯
            "#FFDB4F";  % 26 黄支子
            "#B3ADA0";  % 28 利休白茶
            "#2ADD9C";  % 30 碧绿
            "#2F4F4F";  % 16 深石板灰
        ];
    end

    % -------- 1. 先画边 --------
    if nargin >= 3 && ~isempty(pointE)
        [rows, cols] = find(pointE ~= 0);
        for k = 1:numel(rows)
            i = rows(k);
            j = cols(k);

            xi = data(i, 1);  yi = data(i, 2);
            xj = data(j, 1);  yj = data(j, 2);

            dx = (xj - xi) * 0.85;
            dy = (yj - yi) * 0.85;

            quiver(xi, yi, dx, dy, 0, ...
                'Color', [0.6 0.6 0.6], ...
                'LineWidth', 0.8, ...
                'MaxHeadSize', 0.5, ...
                'AutoScale', 'off', ...
                'HandleVisibility','off');
        end
    end

    % -------- 2. 准备标签信息 --------
    pos_labels = labels(labels > 0);
    uniq_pos   = unique(pos_labels(:));    % 当前图中的簇ID集合
    K          = numel(uniq_pos);

    idx_noise  = (labels <= 0);
    n_noise    = nnz(idx_noise);

    % -------- 3. 按簇绘制 --------
    for i = 1:K
        cid = uniq_pos(i);                 % 当前簇ID（1,2,3,...）
        idx = (labels == cid);

        % 用簇ID决定颜色（保证相同ID用同一颜色）
        color_idx = 1 + mod(cid-1, size(base_cmap,1));
        hex_color = base_cmap(color_idx);
        c = hex2rgb(hex_color);

        scatter(data(idx,1), data(idx,2), 20, ...
                'MarkerFaceColor', c, ...
                'MarkerEdgeColor', 'k', ...
                'DisplayName', sprintf('C%d', cid));
    end

    % -------- 4. 噪声点 --------
    if n_noise > 0
        scatter(data(idx_noise,1), data(idx_noise,2), 20, ...
                'MarkerFaceColor', [0.5 0.5 0.5], ...
                'MarkerEdgeColor', 'none', ...
                'DisplayName', 'Noise');
    end

    axis equal;
    legend('show','Location','bestoutside');
end

function rgb = hex2rgb(hex)
    if isstring(hex)
        hex = char(hex);
    end
    if startsWith(hex, "#")
        hex = hex(2:end);
    end
    rgb = sscanf(hex, '%2x%2x%2x', [1 3]) / 255;
end
