function[A]=demo_Visualization(position_1, position_2, im1, im2, cleanedPoints1, cleanedPoints2, file1, file2,cor1, rmse,H)

A = 1;
matchedPoints1 = cor1(:,1:2);

% 创建保存结果的文件夹 - 使用绝对路径确保位置明确
current_folder = pwd;
output_folder = fullfile(current_folder, '拼接结果');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('已创建保存文件夹: %s\n', output_folder);
else
    fprintf('保存文件夹已存在: %s\n', output_folder);
end

% 获取文件名（不含扩展名）用于保存
[~, name1] = fileparts(file1);
[~, name2] = fileparts(file2);

fprintf('\n=== 开始显示特征点 ===\n');

% 显示第一张图片的特征点
fig_features1 = figure('Name', '图像1特征点', 'NumberTitle', 'off', 'Visible', 'on');
imshow(im1); hold on;
plot(position_1(:,1), position_1(:,2), 'ro', 'MarkerSize', 3, 'LineWidth', 1);
title(['图像1特征点 (' num2str(size(position_1,1)) ' 个)']);
hold off;

% 保存第一张图片的特征点图像
features1_filename = fullfile(output_folder, [name1 '_features.png']);
saveas(fig_features1, features1_filename);
fprintf('✓ 已保存图像1特征点: %s\n', features1_filename);
set(fig_features1, 'Visible', 'on');

% 显示第二张图片的特征点
fig_features2 = figure('Name', '图像2特征点', 'NumberTitle', 'off', 'Visible', 'on');
imshow(im2); hold on;
plot(position_2(:,1), position_2(:,2), 'ro', 'MarkerSize', 3, 'LineWidth', 1);
title(['图像2特征点 (' num2str(size(position_2,1)) ' 个)']);
hold off;

% 保存第二张图片的特征点图像
features2_filename = fullfile(output_folder, [name2 '_features.png']);
saveas(fig_features2, features2_filename);
fprintf('✓ 已保存图像2特征点: %s\n', features2_filename);
set(fig_features2, 'Visible', 'on');



%% ============ 保存特征点匹配图像 ============
fprintf('\n=== 开始保存结果 ===\n');
fig_matched = figure('Name', '特征点匹配连线', 'NumberTitle', 'off', 'Visible', 'off');
showMatchedFeatures(im1, im2, cleanedPoints1, cleanedPoints2, 'montage');
title('匹配点连线');

% 保存特征点匹配图像

matched_points_filename = fullfile(output_folder, [name1 '_' name2 '_matched_points.png']);
saveas(fig_matched, matched_points_filename);
fprintf('✓ 已保存特征点匹配图像: %s\n', matched_points_filename);


% 显示特征点匹配图像
set(fig_matched, 'Visible', 'on');



%%


allpoints = size(matchedPoints1, 1);
CM = size(cleanedPoints1, 1);
CMR = (CM / allpoints) * 100;
fprintf('初始匹配点数: %d\n', allpoints);
fprintf('正确匹配点数 (CM): %d\n', CM);
fprintf('正确匹配率 (CMR): %.2f%%\n', CMR);
disp(['RMSE is : ',num2str(rmse),'pixel']);




%% ============ 图像变换部分 ============
% 设置变换参数
trans_form = 'similarity'; % 变换类型：'similarity', 'affine', 'projective', 'polynomial-2'
out_form = 'union';        % 输出范围：'reference', 'union', 'inter', 'geo'
chg_scale = 1;             % 是否改变尺度
show_flag = 1;             % 是否显示变换结果
overlap_flag = 1;          % 是否生成重叠图像
mosaic_flag = 1;           % 是否生成拼接图像

% 调用图像变换函数
fprintf('开始图像变换...\n');
tic;
[~, ~, I1_rs, I2_rs, I3, I4, ~, ~] = Transformation(im1, im2, cleanedPoints1, cleanedPoints2, trans_form, out_form, chg_scale, show_flag, overlap_flag, mosaic_flag);
t_trans = toc; 
fprintf('已完成图像变换，用时 %.2fs\n', t_trans);

% ============ 保存重叠图像和拼接图像 ============
if overlap_flag && ~isempty(I3)
    % 显示重叠图像
    fig_overlap = figure('Name', '重叠图像', 'NumberTitle', 'off');
    imshow(I3); title('Overlap Form'); drawnow
    
    % 保存重叠图像
    overlap_filename = fullfile(output_folder, [name1 '_' name2 '_overlap.png']);
    imwrite(I3, overlap_filename);
    fprintf('✓ 已保存重叠图像: %s\n', overlap_filename);
    
    close(fig_overlap); % 关闭图形窗口以避免太多窗口
end

if mosaic_flag && ~isempty(I4)
    % 显示拼接图像
    fig_mosaic = figure('Name', '拼接图像', 'NumberTitle', 'off');
    imshow(I4); title('Mosaic Form'); drawnow
    
    % 保存拼接图像
    mosaic_filename = fullfile(output_folder, [name1 '_' name2 '_mosaic.png']);
    imwrite(I4, mosaic_filename);
    fprintf('✓ 已保存拼接图像: %s\n', mosaic_filename);
    
    close(fig_mosaic); % 关闭图形窗口
end



% ============ 保存变换后的单独图像 ============
if ~isempty(I1_rs)
    transformed1_filename = fullfile(output_folder, [name1 '_transformed.png']);
    imwrite(I1_rs, transformed1_filename);
    %fprintf('✓ 已保存变换后的图像1: %s\n', transformed1_filename);
end

if ~isempty(I2_rs)
    transformed2_filename = fullfile(output_folder, [name2 '_transformed.png']);
    imwrite(I2_rs, transformed2_filename);
    %fprintf('✓ 已保存变换后的图像2: %s\n', transformed2_filename);
end

% ============ 保存匹配结果信息到文本文件 ============
info_filename = fullfile(output_folder, [name1 '_' name2 '_matching_info.txt']);
fid = fopen(info_filename, 'w');
if fid ~= -1
    fprintf(fid, '图像拼接匹配结果信息\n');
    fprintf(fid, '=====================\n\n');
    fprintf(fid, '输入图像1: %s\n', file1);
    fprintf(fid, '输入图像2: %s\n\n', file2);
    fprintf(fid,'The number of feature points from the reference image is: %d.\n', size(position_1,1));
    fprintf(fid,'The number of feature points from the image to be registered is: %d.\n', size(position_2,1));
    fprintf(fid, '匹配统计:\n');
    fprintf(fid, '  初始匹配点数: %d\n', allpoints);
    fprintf(fid, '  正确匹配点数 (CM): %d\n', CM);
    fprintf(fid, '  正确匹配率 (CMR): %.2f%%\n', CMR);
    fprintf(fid, '  RMSE: %.4f 像素\n\n', rmse);
    % fprintf(fid, '时间统计:\n');
    % fprintf(fid, '  特征匹配时间: %.4f 秒\n', time);
    % fprintf(fid, '  图像变换时间: %.4f 秒\n', t_trans);
    % fprintf(fid, '  总处理时间: %.4f 秒\n\n', time + t_trans);
    fprintf(fid, '变换矩阵 H:\n');
    for i = 1:size(H,1)
        fprintf(fid, '  [');
        fprintf(fid, '%.6f ', H(i,:));
        fprintf(fid, ']\n');
    end
    fclose(fid);
    fprintf('✓ 已保存匹配信息: %s\n', info_filename);
end


% 尝试打开保存文件夹（Windows系统）
if ispc
    try
        winopen(output_folder);
        fprintf('已自动打开保存文件夹\n');
    catch
        fprintf('请手动打开文件夹: %s\n', output_folder);
    end
else
    fprintf('请手动打开文件夹: %s\n', output_folder);
end

% ============ 图像变换部分结束 ============

% save RES_cofsm.mat RES



end