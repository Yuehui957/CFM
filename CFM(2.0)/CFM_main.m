function CFM_main()


clc;
close all;

warning('off')


%% =========================
% 添加所有算法路径
% =========================


root=pwd;



addpath(fullfile(root,'MI'));

addpath(fullfile(root,'Visualization'));

addpath(fullfile(root,'1ScaleSpaceConstruction'));

addpath(fullfile(root,'2KeypointDetection'));

addpath(fullfile(root,'3FeatureDescriptors'));

addpath(fullfile(root,'4FeatureMatching'));



%% =========================
% 调用你的原始程序
% =========================


demo_CFM_core();



end