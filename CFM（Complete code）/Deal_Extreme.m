function [I_deal, flag] = Deal_Extreme(I, method, ex_size, flag)
% 处理极端值
    if nargin < 2
        method = 1;
    end
    if nargin < 3
        ex_size = 4096;
    end
    if nargin < 4
        flag = 1;
    end
    
    [M, N, B] = size(I);
    if max(M, N) > ex_size && flag
        if M > N
            N = round(N * ex_size / M);
            M = ex_size;
        else
            M = round(M * ex_size / N);
            N = ex_size;
        end
        I_deal = imresize(I, [M, N]);
        flag = 0;
    else
        I_deal = I;
    end
end

function I_vis = Visual(I)
% 可视化图像
    if isinteger(I)
        I = double(I);
    end
    I_vis = I - min(I(:));
    I_vis = I_vis / max(I_vis(:));
    I_vis = uint8(I_vis * 255);
end

function [I1_out, I2_out, mosaic] = Mosaic_Map(I1, I2, grid_size)
% 创建马赛克拼接图
    [M, N, B] = size(I1);
    mosaic = zeros(M, N, B);
    
    for i = 1:grid_size:M
        for j = 1:grid_size:N
            if mod(floor(i/grid_size) + floor(j/grid_size), 2) == 0
                i_end = min(i+grid_size-1, M);
                j_end = min(j+grid_size-1, N);
                mosaic(i:i_end, j:j_end, :) = I1(i:i_end, j:j_end, :);
            else
                i_end = min(i+grid_size-1, M);
                j_end = min(j+grid_size-1, N);
                mosaic(i:i_end, j:j_end, :) = I2(i:i_end, j:j_end, :);
            end
        end
    end
    
    I1_out = I1;
    I2_out = I2;
end