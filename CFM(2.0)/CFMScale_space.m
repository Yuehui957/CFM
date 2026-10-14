function [CoOcurscale_space, K_abs] = CFMScale_space(im, sigma_1, Mmax, ratio, sigma_s)

    if isa(im, 'uint8')
        im = double(im);
    end
    [M, N] = size(im);
    CoOcurscale_space = zeros(M, N, Mmax);
    K_abs = zeros(M, N, Mmax);

    sigma = sigma_1 * ratio.^(0:Mmax-1);

    lambda = 1.0;      % λ = 1
    alpha = 0.01;      % 法向扩散系数
    beta = 1.0;        % 切向扩散系数
    dt = 0.1;          % 时间步长
    n_iter = 5;        % 每个尺度的迭代次数（总时间 T = n_iter * dt）


    I_current = im;

    for s = 1:Mmax
        cur_sigma = sigma(s);
        fprintf('第 %d 层尺度空间构建中（sigma = %.3f）...\n', s, cur_sigma);

        % 1. 计算当前尺度下的高斯曲率
        % 先对图像做高斯平滑
        win_size = 2 * ceil(3 * cur_sigma) + 1;
        g = fspecial('gaussian', [win_size, win_size], cur_sigma);
        I_smooth = imfilter(I_current, g, 'replicate', 'same');
        

        [Ix, Iy, Ixx, Ixy, Iyy] = gaussian_derivatives(I_smooth, cur_sigma);


        % 高斯曲率 K
        E = 1 + Ix.^2;
        F = Ix .* Iy;
        G = 1 + Iy.^2;
        denom = sqrt(1 + Ix.^2 + Iy.^2);
        denom(denom < eps) = eps;
        L = Ixx ./ denom;
        M_cur = Ixy ./ denom;
        N_cur = Iyy ./ denom;
        denom_K = E.*G - F.^2;
        denom_K(denom_K < eps) = eps;
        K = (L.*N_cur - M_cur.^2) ./ denom_K;
        K_abs_cur = abs(K);
        K_abs(:,:,s) = K_abs_cur;
        
        % 2. 各向异性扩散
        
        I_cur = I_current;
        
        for iter = 1:n_iter
 
            [Ix_cur, Iy_cur] = gradient(I_cur);
            grad_mag = sqrt(Ix_cur.^2 + Iy_cur.^2) + eps;

            nx = Ix_cur ./ grad_mag;
            ny = Iy_cur ./ grad_mag;

            tx = -Iy_cur ./ grad_mag;
            ty = Ix_cur ./ grad_mag;
            

            D = 1 ./ (1 + lambda * K_abs_cur);
            


            Dxx = D .* (alpha * nx.^2 + beta * tx.^2);
            Dxy = D .* (alpha * nx.*ny + beta * tx.*ty);
            Dyy = D .* (alpha * ny.^2 + beta * ty.^2);
            


            Fx = Dxx .* Ix_cur + Dxy .* Iy_cur;
            Fy = Dxy .* Ix_cur + Dyy .* Iy_cur;
            

            div_Fx = zeros(size(I_cur));
            div_Fy = zeros(size(I_cur));
           
            div_Fx(:, 2:end-1) = (Fx(:, 3:end) - Fx(:, 1:end-2)) / 2;
            div_Fy(2:end-1, :) = (Fy(3:end, :) - Fy(1:end-2, :)) / 2;
            
            div_Fx(:, 1) = Fx(:, 2) - Fx(:, 1);
            div_Fx(:, end) = Fx(:, end) - Fx(:, end-1);
            div_Fy(1, :) = Fy(2, :) - Fy(1, :);
            div_Fy(end, :) = Fy(end, :) - Fy(end-1, :);
            
            div = div_Fx + div_Fy;
          
            I_cur = I_cur + dt * div;

        end


        CoOcurscale_space(:,:,s) = I_cur;

        I_current = I_cur;
    end
end

function [Ix, Iy, Ixx, Ixy, Iyy] = gaussian_derivatives(I, sigma)
%%   
    h1 = [-1,0,1;-2,0,2;-1,0,1];
    h2 = [-1,-2,-1;0,0,0;1,2,1];
    lp = lowpassfilter([3,3], .45, 20);


    gradient_x_0 = imfilter(I, lp, 'replicate');
    gradient_y_0 = imfilter(I, lp, 'replicate');
    gradient_x_0 = real(gradient_x_0);
    gradient_y_0 = real(gradient_y_0);
    gradient_0 = sqrt(gradient_x_0.^2 + gradient_y_0.^2);

    gradient_x_1 = imfilter(gradient_0, h1, 'replicate');
    gradient_y_1 = imfilter(gradient_0, h2, 'replicate');
    I0 = sqrt(gradient_x_1.^2 + gradient_y_1.^2);


%%
    Ix = imfilter(I0, h1, 'replicate');
    Iy = imfilter(I0, h2, 'replicate');
    I1 = sqrt(Ix.^2 + Iy.^2);

    Ixx = imfilter(I1, h1, 'replicate');
    Iyy = imfilter(I1, h2, 'replicate');
    Ixy = sqrt(Ixx.^2 + Iyy.^2);
    % Ixy = imfilter(Ix, h2, 'replicate');


end