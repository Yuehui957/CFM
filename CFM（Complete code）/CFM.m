function[H,rmse,cor1,cor2,position_1,position_2]=CFM(CoF1,CoF2,image1, image2,d_SH,sigma_1,ratio,first_layer,change_form,K_abs1,K_abs2)
    t1=clock;
    tic;
    

    [CFM_1, gradient_1, angle_1] = Shi_tomasi_scale_improve(CoF1, sigma_1, ratio, K_abs1);
    [CFM_2, gradient_2, angle_2] = Shi_tomasi_scale_improve(CoF2, sigma_1, ratio, K_abs2);
    disp(['The CFM scale space spend time is：',num2str(toc),'S']);
    tic;
  
    % 特征点提取
    [position_1] = find_scale_extreme(CFM_1, d_SH, sigma_1, ratio, gradient_1, angle_1, first_layer,image1);
    [position_2] = find_scale_extreme(CFM_2, d_SH, sigma_1, ratio, gradient_2, angle_2, first_layer,image2);
    
    fprintf('The number of feature points from the reference image is: %d.\n', size(position_1,1));
    fprintf('The number of feature points from the image to be registered is: %d.\n', size(position_2,1));
    
    disp(['The time taken to extract feature points is：',num2str(toc),'S']);
    tic;
    
                                             
    [descriptors_1, locs_1] = GLOH_descriptors_improve(gradient_1, angle_1, position_1,K_abs1);
    [descriptors_2, locs_2] = GLOH_descriptors_improve( gradient_2, angle_2, position_2,K_abs2);


    disp(['The time taken to compute the descriptor is：',num2str(toc),'S']);
    tic;
  
    [solution, rmse, cor1, cor2] = match(descriptors_2, locs_2, descriptors_1, locs_1, change_form);
    H = inv(solution);
    
    t2 = clock;
    disp(['The feature points extraction and matching time spent are：',num2str(etime(t2,t1)),'S']);




end