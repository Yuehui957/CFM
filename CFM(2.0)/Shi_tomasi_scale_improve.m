function [MSC_R_function, gradient, angle] = Shi_tomasi_scale_improve(CFMscale_space, sigma, ratio, precomputed_K)


    if nargin < 4
        precomputed_K = [];
    end

    [M, N, P] = size(CFMscale_space);
    MSC_R_function = zeros(M, N, P);
    gradient = zeros(M, N, P);
    angle = zeros(M, N, P);
    
    % 梯度计算核
    h1 = [-1,0,1;-2,0,2;-1,0,1];
    h2 = [-1,-2,-1;0,0,0;1,2,1];
    lp = lowpassfilter([3,3], .45, 20);
    
    K_weight_all = cell(1, P);
    global_max_K = 0;
    
    for i = 1:P
        I = CFMscale_space(:,:,i);
        
        % ----- 梯度与角度计算-----
        gradient_x_0 = imfilter(I, lp, 'replicate');
        gradient_y_0 = imfilter(I, lp, 'replicate');
        gradient_x_0 = real(gradient_x_0);
        gradient_y_0 = real(gradient_y_0);
        gradient_0 = sqrt(gradient_x_0.^2 + gradient_y_0.^2);
        
        gradient_x_1 = imfilter(gradient_0, h1, 'replicate');
        gradient_y_1 = imfilter(gradient_0, h2, 'replicate');
        gradient_1 = sqrt(gradient_x_1.^2 + gradient_y_1.^2);
        
        gradient_x_2 = imfilter(gradient_1, h1, 'replicate');
        gradient_y_2 = imfilter(gradient_1, h2, 'replicate');
        gradient_2 = sqrt(gradient_x_2.^2 + gradient_y_2.^2);
        angle_2 = atan2(gradient_y_2, gradient_x_2);
        angle_2 = angle_2 * 180 / pi;
        angle_2(angle_2 < 0) = angle_2(angle_2 < 0) + 360;
        

        % ----- 结构张量计算-----
        cur_scale = sigma * ratio^(i-1);
        Csh_11 = cur_scale^2 * gradient_x_1 .* gradient_x_1;
        Csh_12 = cur_scale^2 * gradient_x_1 .* gradient_y_1;
        Csh_22 = cur_scale^2 * gradient_y_1 .* gradient_y_1;
        
        gaussian_sigma = cur_scale * sqrt(2);
        temp = round(2 * gaussian_sigma);
        gaussian_width = 2 * temp + 1;
        W = fspecial('gaussian', [gaussian_width, gaussian_width], gaussian_sigma);
        [a,b] = meshgrid(1:gaussian_width, 1:gaussian_width);
        index = ((a-temp-1).^2 + (b-temp-1).^2) > temp^2;
        W(index) = 0;
        
        Csh_11 = imfilter(Csh_11, W, 'replicate');
        Csh_12 = imfilter(Csh_12, W, 'replicate');
        Csh_22 = imfilter(Csh_22, W, 'replicate');
        
        T = Csh_11 + Csh_22;
        D = Csh_11 .* Csh_22 - Csh_12.^2;
        R_original = T/2 - sqrt(T.^2/4 - D);  % 原始结构张量响应
        
        % ----- 曲率权重获取 -----
        if ~isempty(precomputed_K) && size(precomputed_K,3) == P
           
            K_weight = precomputed_K(:,:,i);
        else
            
            I_double = double(I);
            cur_sigma = sigma * ratio^(i-1);
            [Ix, Iy, Ixx, Ixy, Iyy] = gaussian_derivatives(I_double, cur_sigma);
            
            E = 1 + Ix.^2;
            F = Ix .* Iy;
            G = 1 + Iy.^2;
            
            denom = sqrt(1 + Ix.^2 + Iy.^2);
            denom(denom < eps) = eps;
            L = Ixx ./ denom;
            M_cur = Ixy ./ denom;
            N_cur = Iyy ./ denom;
            
            denom_K = E.*G - F.^2;
            K = zeros(size(denom_K));
            valid = denom_K > 1e-8;
            K(valid) = (L(valid).*N_cur(valid) - M_cur(valid).^2) ./ denom_K(valid);
            K_weight = abs(K);
        end
        
        K_weight_all{i} = K_weight;
        global_max_K = max(global_max_K, max(K_weight(:)));
        MSC_R_function(:,:,i) = 100 .* K_weight;
        % MSC_R_function(:,:,i) = R_original;

        gradient(:,:,i) = gradient_2;
        angle(:,:,i) = angle_2;  

    end
    
    fprintf('曲率权重最大值: %.4f\n', global_max_K);
    fprintf('最终响应最大值: %.2f\n', max(MSC_R_function(:)));

    
end




% -------------------------------------------------------------------------
function [Ix, Iy, Ixx, Ixy, Iyy] = gaussian_derivatives(I, sigma)
    r = ceil(3*sigma);
    if r < 1, r = 1; end
    x = -r:r;
    y = x';
    
    G   = exp(-(x.^2)/(2*sigma^2));
    Gx  = -x/sigma^2 .* G;
    Gxx = (x.^2/sigma^4 - 1/sigma^2) .* G;
    
    G_y = exp(-(y.^2)/(2*sigma^2));
    Gy  = -y/sigma^2 .* G_y;
    Gyy = (y.^2/sigma^4 - 1/sigma^2) .* G_y;
    
    Ix   = conv2(conv2(I, Gx', 'same'), G_y, 'same');
    Iy   = conv2(conv2(I, G',  'same'), Gy,  'same');
    Ixx  = conv2(conv2(I, Gxx','same'), G_y, 'same');
    Iyy  = conv2(conv2(I, G',  'same'), Gyy, 'same');
    Ixy  = conv2(conv2(I, Gx', 'same'), Gy,  'same');
end