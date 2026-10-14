function [S, RMO, delat_x, delat_y, scale_hist, scale_x, angle_hist, ...
          angle_x, x_hist, x_x, y_hist, y_x] = XY_position_offset(cor1, cor2, unit_scale, unit_angle, unit_x, unit_y)
    
    % 初始化所有输出变量，防止未定义错误
    S = 1;
    RMO = 0;
    delat_x = 0;
    delat_y = 0;
    scale_hist = [];
    scale_x = [];
    angle_hist = [];
    angle_x = [];
    x_hist = [];
    x_x = [];
    y_hist = [];
    y_x = [];
    
    % 安全索引辅助函数
    function result = safe_index(array, indices)
        if isempty(indices) || isempty(array)
            result = [];
        else
            try
                result = array(indices);
            catch
                result = [];
            end
        end
    end
    
    % 检查输入有效性
    if isempty(cor1) || isempty(cor2) || size(cor1,1) < 2 || size(cor2,1) < 2
        return;
    end
    
    smooth = 'true';
    
    %% 1. 处理尺度比
    ratio_cor_scale = cor1(:,1) ./ cor2(:,1);
    ratio_cor_scale = ratio_cor_scale';
    min_scale = min(ratio_cor_scale);
    max_scale = max(ratio_cor_scale);
    
    scale_bin = round((max_scale - min_scale) / unit_scale);
    if scale_bin < 0
        scale_bin = 0;
    end
    
    scale_hist = zeros(1, scale_bin + 1);
    ratio_cor_scale_1 = round((ratio_cor_scale - min_scale) / unit_scale);
    [~, num_scale] = size(ratio_cor_scale_1);
    
    for i = 1:num_scale
        idx = ratio_cor_scale_1(1,i) + 1;
        if idx >= 1 && idx <= length(scale_hist)
            scale_hist(idx) = scale_hist(idx) + 1;
        end
    end
    
    if strcmp(smooth, 'true')
        n = scale_bin + 1;
        if n >= 5
            hist = zeros(1, n);
            if n > 1
                hist(1) = (scale_hist(n-1) + scale_hist(3))/16 + 4*(scale_hist(n) + scale_hist(2))/16 + scale_hist(1)*6/16;
            end
            if n > 2
                hist(2) = (scale_hist(n) + scale_hist(4))/16 + 4*(scale_hist(1) + scale_hist(3))/16 + scale_hist(2)*6/16;
            end
            if n > 4
                hist(3:n-2) = (scale_hist(1:n-4) + scale_hist(5:n))/16 + 4*(scale_hist(2:n-3) + scale_hist(4:n-1))/16 + scale_hist(3:n-2)*6/16;
            end
            if n > 2
                hist(n-1) = (scale_hist(n-3) + scale_hist(1))/16 + 4*(scale_hist(n-2) + scale_hist(n))/16 + scale_hist(n-1)*6/16;
            end
            if n > 1
                hist(n) = (scale_hist(n-2) + scale_hist(2))/16 + 4*(scale_hist(n-1) + scale_hist(1))/16 + scale_hist(n)*6/16;
            end
            scale_hist = hist;
        end
    end
    
    x_length = size(scale_hist, 2);
    if x_length > 0
        scale_x = zeros(1, x_length);
        for i = 1:x_length
            scale_x(1,i) = (i-1) * unit_scale + min_scale;
        end
    end
    
    % 应用过滤器
    if ~isempty(scale_x)
        index = find(scale_x < 5);
        scale_x = safe_index(scale_x, index);
        scale_hist = safe_index(scale_hist, index);
    end
    
    % 找到最大值的索引
    if ~isempty(scale_hist)
        [~, index_scale] = max(scale_hist);
    else
        index_scale = 1;
    end
    
    if isempty(index_scale)
        index_scale = 1;
    end
    
    n = length(scale_hist);
    
    if strcmp(smooth, 'true') && n > 0
        k1 = index_scale - 1;
        k2 = index_scale + 1;
        
        if k1 < 1
            k1 = n;
        end
        if k2 > n
            k2 = 1;
        end
        
        if n > 1
            denominator = scale_hist(k1) + scale_hist(k2) - 2 * scale_hist(index_scale) + 10^-5;
            if abs(denominator) > eps
                index_scale = index_scale + 0.5 * (scale_hist(k1) - scale_hist(k2)) / denominator;
            end
        end
    end
    
    if ~isempty(scale_x) && ~isempty(index_scale)
        S = (index_scale - 1) * unit_scale + scale_x(1);
    else
        S = 1;
    end
    S = mean(S);
    
    %% 2. 处理角度差
    diff_cor_angle = cor1(:,2) - cor2(:,2);
    diff_cor_angle = diff_cor_angle';
    min_diff_cor_angle = min(diff_cor_angle);
    max_diff_cor_angle = max(diff_cor_angle);
    
    angle_bin = round((max_diff_cor_angle - min_diff_cor_angle) / unit_angle);
    if angle_bin < 0
        angle_bin = 0;
    end
    
    angle_hist = zeros(1, angle_bin + 1);
    diff_cor_angle_1 = round((diff_cor_angle - min_diff_cor_angle) / unit_angle);
    [~, num] = size(diff_cor_angle_1);
    
    for i = 1:num
        idx = diff_cor_angle_1(1,i) + 1;
        if idx >= 1 && idx <= length(angle_hist)
            angle_hist(idx) = angle_hist(idx) + 1;
        end
    end
    
    if strcmp(smooth, 'true')
        n = angle_bin + 1;
        if n >= 5
            hist = zeros(1, n);
            if n > 1
                hist(1) = (angle_hist(n-1) + angle_hist(3))/16 + 4*(angle_hist(n) + angle_hist(2))/16 + angle_hist(1)*6/16;
            end
            if n > 2
                hist(2) = (angle_hist(n) + angle_hist(4))/16 + 4*(angle_hist(1) + angle_hist(3))/16 + angle_hist(2)*6/16;
            end
            if n > 4
                hist(3:n-2) = (angle_hist(1:n-4) + angle_hist(5:n))/16 + 4*(angle_hist(2:n-3) + angle_hist(4:n-1))/16 + angle_hist(3:n-2)*6/16;
            end
            if n > 2
                hist(n-1) = (angle_hist(n-3) + angle_hist(1))/16 + 4*(angle_hist(n-2) + angle_hist(n))/16 + angle_hist(n-1)*6/16;
            end
            if n > 1
                hist(n) = (angle_hist(n-2) + angle_hist(2))/16 + 4*(angle_hist(n-1) + angle_hist(1))/16 + angle_hist(n)*6/16;
            end
            angle_hist = hist;
        end
    end
    
    y_length = size(angle_hist, 2);
    if y_length > 0
        angle_x = zeros(1, y_length);
        for i = 1:y_length
            angle_x(1,i) = (i-1) * unit_angle + min_diff_cor_angle;
        end
    end
    
    % 找到最大值的索引
    if ~isempty(angle_hist)
        [~, index_angle] = max(angle_hist);
    else
        index_angle = 1;
    end
    
    if isempty(index_angle)
        index_angle = 1;
    end
    
    n = length(angle_hist);
    
    if strcmp(smooth, 'true') && n > 0
        k1 = index_angle - 1;
        k2 = index_angle + 1;
        
        if k1 < 1
            k1 = n;
        end
        if k2 > n
            k2 = 1;
        end
        
        if n > 1
            denominator = angle_hist(k1) + angle_hist(k2) - 2 * angle_hist(index_angle) + 10^-5;
            if abs(denominator) > eps
                index_angle = index_angle + 0.5 * (angle_hist(k1) - angle_hist(k2)) / denominator;
            end
        end
    end
    
    if ~isempty(angle_x) && ~isempty(index_angle)
        RMO = (index_angle - 1) * unit_angle + angle_x(1);
    else
        RMO = 0;
    end
    RMO = mean(RMO);
    
    %% 3. 处理X方向偏移
    cor1_xy = cor1(:, [3,4]);
    cor2_xy = cor2(:, [3,4]);
    
    temp_S = 1 / S;
    temp_RMO = -RMO;
    temp_T = [temp_S * cos(temp_RMO/180*pi), -temp_S * sin(temp_RMO/180*pi);
              temp_S * sin(temp_RMO/180*pi), temp_S * cos(temp_RMO/180*pi)];
    
    diff_xy = temp_T * cor1_xy' - cor2_xy';
    diff_x = diff_xy(1,:);
    
    min_diff_x = min(diff_x);
    max_diff_x = max(diff_x);
    
    x_bin = round((max_diff_x - min_diff_x) / unit_x);
    if x_bin < 0
        x_bin = 0;
    end
    
    x_hist = zeros(1, x_bin + 1);
    diff_x_1 = round((diff_x - min_diff_x) / unit_x);
    [~, num] = size(diff_x_1);
    
    for i = 1:num
        idx = diff_x_1(1,i) + 1;
        if idx >= 1 && idx <= length(x_hist)
            x_hist(idx) = x_hist(idx) + 1;
        end
    end
    
    if strcmp(smooth, 'true')
        n = x_bin + 1;
        if n >= 5
            hist = zeros(1, n);
            if n > 1
                hist(1) = (x_hist(n-1) + x_hist(3))/16 + 4*(x_hist(n) + x_hist(2))/16 + x_hist(1)*6/16;
            end
            if n > 2
                hist(2) = (x_hist(n) + x_hist(4))/16 + 4*(x_hist(1) + x_hist(3))/16 + x_hist(2)*6/16;
            end
            if n > 4
                hist(3:n-2) = (x_hist(1:n-4) + x_hist(5:n))/16 + 4*(x_hist(2:n-3) + x_hist(4:n-1))/16 + x_hist(3:n-2)*6/16;
            end
            if n > 2
                hist(n-1) = (x_hist(n-3) + x_hist(1))/16 + 4*(x_hist(n-2) + x_hist(n))/16 + x_hist(n-1)*6/16;
            end
            if n > 1
                hist(n) = (x_hist(n-2) + x_hist(2))/16 + 4*(x_hist(n-1) + x_hist(1))/16 + x_hist(n)*6/16;
            end
            x_hist = hist;
        end
    end
    
    y_length = size(x_hist, 2);
    if y_length > 0
        x_x = zeros(1, y_length);
        for i = 1:y_length
            x_x(1,i) = (i-1) * unit_x + min_diff_x;
        end
    end
    
    % 应用过滤器
    if ~isempty(x_x)
        index = find(x_x > -1000 & x_x < 1000);
        x_x = safe_index(x_x, index);
        x_hist = safe_index(x_hist, index);
    end
    
    % 找到最大值的索引
    if ~isempty(x_hist)
        [~, index_x] = max(x_hist);
    else
        index_x = 1;
    end
    
    if isempty(index_x)
        index_x = 1;
    end
    
    n = length(x_hist);
    
    if strcmp(smooth, 'true') && n > 0
        k1 = index_x - 1;
        k2 = index_x + 1;
        
        if k1 < 1
            k1 = n;
        end
        if k2 > n
            k2 = 1;
        end
        
        if n > 1
            denominator = x_hist(k1) + x_hist(k2) - 2 * x_hist(index_x) + 10^-5;
            if abs(denominator) > eps
                index_x = index_x + 0.5 * (x_hist(k1) - x_hist(k2)) / denominator;
            end
        end
    end
    
    if ~isempty(x_x) && ~isempty(index_x)
        delat_x = (index_x - 1) * unit_x + x_x(1);
    else
        delat_x = 0;
    end
    delat_x = mean(delat_x);
    
    %% 4. 处理Y方向偏移
    diff_y = diff_xy(2,:);
    min_diff_y = min(diff_y);
    max_diff_y = max(diff_y);
    
    y_bin = round((max_diff_y - min_diff_y) / unit_y);
    if y_bin < 0
        y_bin = 0;
    end
    
    y_hist = zeros(1, y_bin + 1);
    diff_y_1 = round((diff_y - min_diff_y) / unit_y);
    [~, num] = size(diff_y_1);
    
    for i = 1:num
        idx = diff_y_1(1,i) + 1;
        if idx >= 1 && idx <= length(y_hist)
            y_hist(idx) = y_hist(idx) + 1;
        end
    end
    
    if strcmp(smooth, 'true')
        n = y_bin + 1;
        if n >= 5
            hist = zeros(1, n);
            if n > 1
                hist(1) = (y_hist(n-1) + y_hist(3))/16 + 4*(y_hist(n) + y_hist(2))/16 + y_hist(1)*6/16;
            end
            if n > 2
                hist(2) = (y_hist(n) + y_hist(4))/16 + 4*(y_hist(1) + y_hist(3))/16 + y_hist(2)*6/16;
            end
            if n > 4
                hist(3:n-2) = (y_hist(1:n-4) + y_hist(5:n))/16 + 4*(y_hist(2:n-3) + y_hist(4:n-1))/16 + y_hist(3:n-2)*6/16;
            end
            if n > 2
                hist(n-1) = (y_hist(n-3) + y_hist(1))/16 + 4*(y_hist(n-2) + y_hist(n))/16 + y_hist(n-1)*6/16;
            end
            if n > 1
                hist(n) = (y_hist(n-2) + y_hist(2))/16 + 4*(y_hist(n-1) + y_hist(1))/16 + y_hist(n)*6/16;
            end
            y_hist = hist;
        end
    end
    
    y_length = size(y_hist, 2);
    if y_length > 0
        y_x = zeros(1, y_length);
        for i = 1:y_length
            y_x(1,i) = (i-1) * unit_y + min_diff_y;
        end
    end
    
    % 应用过滤器
    if ~isempty(y_x)
        index = find(y_x > -1000 & y_x < 1000);
        y_x = safe_index(y_x, index);
        y_hist = safe_index(y_hist, index);
    end
    
    % 找到最大值的索引
    if ~isempty(y_hist)
        [~, index_y] = max(y_hist);
    else
        index_y = 1;
    end
    
    if isempty(index_y)
        index_y = 1;
    end
    
    n = length(y_hist);
    
    if strcmp(smooth, 'true') && n > 0
        k1 = index_y - 1;
        k2 = index_y + 1;
        
        if k1 < 1
            k1 = n;
        end
        if k2 > n
            k2 = 1;
        end
        
        if n > 1
            denominator = y_hist(k1) + y_hist(k2) - 2 * y_hist(index_y) + 10^-5;
            if abs(denominator) > eps
                index_y = index_y + 0.5 * (y_hist(k1) - y_hist(k2)) / denominator;
            end
        end
    end
    
    if ~isempty(y_x) && ~isempty(index_y)
        delat_y = (index_y - 1) * unit_y + y_x(1);
    else
        delat_y = 0;
    end
    delat_y = mean(delat_y);
    
end