clc;clear;close all; warning('off')



RES=[];



[file1, path1] = uigetfile({'*.jpg;*.png;*.bmp;*.tif','Image Files'}, '选择第一张图像');
if isequal(file1,0)
    error('未选择图像');
end
str1 = fullfile(path1, file1);

[file2, path2] = uigetfile({'*.jpg;*.png;*.bmp;*.tif','Image Files'}, '选择第二张图像');
if isequal(file2,0)
    error('未选择图像');
end
str2 = fullfile(path2, file2);



%% 
% im1=imread("pair11-1.jpg");
% im2=imread("pair11-2.jpg");

im1 = uint8(imread(str1));
im2 = uint8(imread(str2));

[~,~,num1]=size(im1);
[~,~,num2]=size(im2);
if(num1==3)
    image_11=rgb2gray(im1);
else
    image_11=im1;
end
if(num2==3)
    image_22=rgb2gray(im2);
else
    image_22=im2;
end

% [gtfile, gtpath] = uigetfile('*.txt', '选择真值文件');
% if isequal(gtfile,0)
%      error('未选择真值文件');
%  end
%  gtstr = fullfile(gtpath, gtfile);
%   gt=load(gtstr);


t1=clock;
[H,rmse,cor2,cor1,position_1,position_2] = CFM_Match(image_11,image_22);
t2=clock;
time=etime(t2,t1);

matchedPoints1 = cor1(:,1:2);
matchedPoints2 = cor2(:,1:2);

% H=[gt;0 0 1];
Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
inliersIndex=E<3;
cleanedPoints1 = matchedPoints1(inliersIndex, :);
cleanedPoints2 = matchedPoints2(inliersIndex, :);
[cleanedPoints2,IA] = unique(cleanedPoints2,'rows');
cleanedPoints1 = cleanedPoints1(IA,:);
cleanedPoints=[cleanedPoints1 cleanedPoints2];
cleanedPoints = double(cleanedPoints);
Y_=H*[cleanedPoints(:,1:2)';ones(1,size(cleanedPoints,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-cleanedPoints(:,3:4)').^2));
if length(E)<10
    rmse1 = 20;
else
    rmse1 = sqrt(sum(E.^2)/size(E,2));
end
rmse = max(rmse,rmse1);
length(E)
timeres = double([time rmse size(cleanedPoints,1)]);
RES = [RES;timeres];
disp(['RMSE is : ',num2str(rmse),'pixel']);
disp(['RMSE1 is : ',num2str(rmse1),'pixel']);

[A] = demo_Visualization(position_1, position_2, im1, im2, cleanedPoints1, cleanedPoints2, file1, file2, cor1,rmse,H);

