function[descriptor] = improved_grid_descriptor_1(gradient, angle, x, y, scale, main_angle, d, n,K_abs,key_point_array)
    cos_t = cos(-main_angle/180*pi);
    sin_t = sin(-main_angle/180*pi);
    [M,N] = size(gradient);
    radius = round(min(12*scale, min(M,N)/3));
    
    radius_x_left = x - radius;
    radius_x_right = x + radius;
    radius_y_up = y - radius;
    radius_y_down = y + radius;
    if(radius_x_left <= 0)
        radius_x_left = 1;
    end
    if(radius_x_right > N)
        radius_x_right = N;
    end
    if(radius_y_up <= 0)
        radius_y_up = 1;
    end
    if(radius_y_down > M)
        radius_y_down = M;
    end

    K_abs = abs(K_abs);
    % 提取局部梯度幅值和原始梯度方向
    sub_gradient = gradient(radius_y_up:radius_y_down, radius_x_left:radius_x_right);
    sub_angle_raw = angle(radius_y_up:radius_y_down, radius_x_left:radius_x_right); 
    sub_K = K_abs(radius_y_up:radius_y_down, radius_x_left:radius_x_right);

    maxK_local = max(sub_K(:));
    if maxK_local > 0
        sub_K_norm = sub_K / maxK_local;
    else
        sub_K_norm = zeros(size(sub_K));
    end

    % ==== 初始化用于补偿的新梯度矩阵和新角度矩阵 ====
    new_gradient = zeros(size(sub_gradient));
    new_angle_raw = zeros(size(sub_angle_raw)); % 补偿后的角度
    
    half_win = 4; % 9x9窗口半径
    [row, col] = size(sub_gradient);

    % if M > 0
    % idx = 1; % 更换为希望展示的实际关键点序号
    % CFM_Fig4(gradient, K_abs, key_point_array(idx,:), ...
    %     fullfile(pwd, 'CFM_Fig4_output'), 'as_provided');
    % end

    for i = 1:row
        for j = 1:col
        % --- 第一部分：计算新的梯度幅值---
            grad_diff_sum = 0;
            neighbor_count = 0;
            center_grad = sub_gradient(i,j);
            for di = -half_win:half_win
                for dj = -half_win:half_win
                    if di==0 && dj==0
                        continue;
                    end
                    ni = i + di;
                    nj = j + dj;
                    if ni>=1 && ni<=row && nj>=1 && nj<=col
                        neighbor_grad = sub_gradient(ni, nj);
                        grad_diff_sum = grad_diff_sum + abs(neighbor_grad - center_grad);
                        neighbor_count = neighbor_count + 1;
                    end
                end
            end
            if neighbor_count > 0
                grad_avg_diff = grad_diff_sum / neighbor_count;   % 平均梯度差异
              
                new_gradient(i,j) = grad_avg_diff ;%* sqrt(sqrt(sub_angle_k(i,j)));
            else
                new_gradient(i,j) = sub_gradient(i,j);
            end
%%            
            % --- 第二部分：基于方向统计学的角度补偿 ---
            center_angle_raw = sub_angle_raw(i, j);
            sum_cos = 0;
            sum_sin = 0;
            neighbor_count_ang = 0;
            
            % 步骤1 & 2: 计算9x9邻域内原始角度的合成向量
            for di = -half_win:half_win
                for dj = -half_win:half_win
                    ni = i + di;
                    nj = j + dj;
                    if ni>=1 && ni<=row && nj>=1 && nj<=col
                        theta_rad = sub_angle_raw(ni, nj) * pi / 180; % 转换为弧度
                        sum_cos = sum_cos + cos(theta_rad);
                        sum_sin = sum_sin + sin(theta_rad);
                        neighbor_count_ang = neighbor_count_ang + 1;
                    end
                end
            end
            
            % 步骤3: 计算平滑后的角度 (邻域平均方向)
            smoothed_angle_rad = atan2(sum_sin, sum_cos);
            % 转换到 [0, 2*pi) 然后到 [0, 360) 度
            if smoothed_angle_rad < 0
                smoothed_angle_rad = smoothed_angle_rad + 2*pi;
            end
            theta_smooth = smoothed_angle_rad * 180 / pi;
            
            % 步骤4: 计算局部方向一致性系数 R
            R = sqrt(sum_cos^2 + sum_sin^2) / neighbor_count_ang;
            
            % 步骤5: 基于一致性系数R，在向量空间进行加权补偿
            % 权重 w = R^2，使得在方向一致的区域更信任平滑结果
            w = R^2;
            
            % 将中心原始角度和平滑角度转换为向量
            v_center = [cos(center_angle_raw*pi/180), sin(center_angle_raw*pi/180)];
            v_smooth = [cos(theta_smooth*pi/180), sin(theta_smooth*pi/180)];
            
            % 计算加权后的最终向量
            v_final = w * v_smooth + (1-w) * v_center;
            
            % 从最终向量解码出补偿后的角度
            theta_comp_rad = atan2(v_final(2), v_final(1));
            if theta_comp_rad < 0
                theta_comp_rad = theta_comp_rad + 2*pi;
            end
            new_angle_raw(i, j) = theta_comp_rad * 180 / pi; % 存储补偿后的原始角度
        end
    end
    
    %% ==== (原)后续处理：使用补偿后的角度进行主方向对齐和量化 ====
    % 注意：这里使用 new_angle_raw 替代原始的 sub_angle_raw
    sub_angle_for_hist = new_angle_raw; 
    sub_angle_for_hist = round((sub_angle_for_hist - main_angle) * n/360);
    sub_angle_for_hist(sub_angle_for_hist <= 0) = sub_angle_for_hist(sub_angle_for_hist <= 0) + n;
    sub_angle_for_hist(sub_angle_for_hist == 0) = n;
    
    % ... [坐标旋转、对数极坐标计算等代码完全保持不变] ...
    X = -(x-radius_x_left):1:(radius_x_right-x);
    Y = -(y-radius_y_up):1:(radius_y_down-y);
    [XX,YY] = meshgrid(X,Y);
    c_rot = XX*cos_t - YY*sin_t;
    r_rot = XX*sin_t + YY*cos_t;
    log_angle = atan2(r_rot, c_rot);
    log_angle = log_angle/pi*180;
    log_angle(log_angle<0) = log_angle(log_angle<0)+360;
    log_amplitude = log2(sqrt(c_rot.^2+r_rot.^2));
    log_angle = round(log_angle*d/360);
    log_angle(log_angle<=0) = log_angle(log_angle<=0)+d;
    log_angle(log_angle>d) = log_angle(log_angle>d)-d;
    r1 = log2(radius*0.25);
    r2 = log2(radius*0.73);
    log_amplitude(log_amplitude<=r1) = 1;
    log_amplitude(log_amplitude>r1&log_amplitude<=r2) = 2;
    log_amplitude(log_amplitude>r2) = 3;
    
    temp_hist = zeros(1, (2*d+1)*n);
    
    % ==== (原)构建直方图：使用 new_gradient 和 sub_angle_for_hist , 后采用曲率====
    for i = 1:row
        for j = 1:col
            angle_bin = log_angle(i,j);
            amplitude_bin = log_amplitude(i,j);
            bin_vertical = sub_angle_for_hist(i,j);
            
            % 原始梯度幅值
            grad_mag = new_gradient(i,j);
            % 曲率权重
            alpha = 0.2;
            k_weight = alpha + (1-alpha) *sub_K_norm(i,j);   % 范围 [0,1]
            
            % 最终权重 = 梯度幅值 * 曲率权重
            Mag = grad_mag * k_weight;
            
            % 直方图累加
            if(amplitude_bin == 1)
                temp_hist(bin_vertical) = temp_hist(bin_vertical) + Mag;
            else
                temp_hist(((amplitude_bin-2)*d+angle_bin-1)*n+bin_vertical+n) = ...
                    temp_hist(((amplitude_bin-2)*d+angle_bin-1)*n+bin_vertical+n) + Mag;
            end
        end
    end
    %%
    temp_hist = temp_hist/sqrt(temp_hist*temp_hist');
    temp_hist(temp_hist>0.2) = 0.2;
    temp_hist = temp_hist/sqrt(temp_hist*temp_hist');
    descriptor = temp_hist;
end