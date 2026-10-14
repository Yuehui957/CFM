function [H,rmse,cor2,cor1,position_1,position_2] = CFM_Match(image_11,image_22)

sigma_s=5;                 % The main parameters of co-occurrence filtering: filter initial window size [3, 5, 10]
sigma_1=1.6;              % Scale of the first layer: 1.6
ratio=2^(1/3);             % scale ratio
Mmax=4;                    % the number of layers in the scale space
first_layer=1;               % extreme point detection start layer number
d_SH=500;                % Feature point extraction threshold size [500, 1000, 1300, 1500, 2000]
change_form='Similarity'; % 'Similarity'£¬'Affine';

fprintf('Reference in image processing:... \n')
[CoOcurscale_space_1,K_abs1]=CFMScale_space(image_11,sigma_1,Mmax,ratio,sigma_s);
fprintf('Image to be registered processing:... \n')
[CoOcurscale_space_2,K_abs2]=CFMScale_space(image_22,sigma_1,Mmax,ratio,sigma_s);

% for i =1:size(CoOcurscale_space_1,3)
%      figure;
%      imshow(CoOcurscale_space_1(:,:,i),[]);
%      title(['³ß¶È¿Õ¼ä-µÚ' num2str(i) '²ã']);
% end

                        
[H, rmse, cor2, cor1,position_1,position_2] = CFM(CoOcurscale_space_1, CoOcurscale_space_2,image_11, image_22, d_SH, sigma_1, ratio, first_layer, change_form,K_abs1,K_abs2);