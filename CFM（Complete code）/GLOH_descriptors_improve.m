function[descriptors,locs]=GLOH_descriptors_improve(gradient,angle,key_point_array,K_abs)
    LOG_POLAR_WIDTH=9;
    LOG_POLAR__HIST_BINS=8;

    M=size(key_point_array,1);
    d=LOG_POLAR_WIDTH;
    n=LOG_POLAR__HIST_BINS;
    descriptors=zeros(M,152);
    locs=key_point_array;
    for i=1:1:M
        x=key_point_array(i,1);
        y=key_point_array(i,2);
        scale=key_point_array(i,3);
        layer=key_point_array(i,4);
        main_angle=key_point_array(i,5);

        current_gradient=gradient(:,:,layer);
        current_angle=angle(:,:,layer);

         descriptors(i,:)=improved_grid_descriptor_1(current_gradient,current_angle,x,y,scale,main_angle,d,n,K_abs);
        
    end
