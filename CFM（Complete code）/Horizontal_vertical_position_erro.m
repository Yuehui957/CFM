function[end_solution,point_1,point_2,rmse]=...
    Horizontal_vertical_position_erro(solution,loc1,loc2,cor1,cor2,des1,des2,Td,change_form,distance)

    loc1=loc1(:,[1,2,5,6,3,4]);
    loc2=loc2(:,[1,2,5,6,3,4]);

    cor1=cor1(:,[3,5,1,2]);
    cor2=cor2(:,[3,5,1,2]);

    unit_scale=0.075;
    unit_angle=9;
    unit_x=7.5;
    unit_y=7.5;
    [S,RMO,delat_x,delat_y,~,~,~,~,~,~,~,~]=XY_position_offset(cor1,cor2,unit_scale,unit_angle,unit_x,unit_y);

    M_des1=size(des1,1);
    M_des2=size(des2,1);

    K1=size(S,2);
    K2=size(RMO,2);
    match=zeros(M_des1,2);
    for i=1:1:K1
        for j=1:1:K2
            prev_match=match;
            for k=1:1:M_des1

                temp1_xy=loc1(k,[1,2]);
                temp1_xy=repmat(temp1_xy,M_des2,1);
                temp1_xy=temp1_xy';
                temp1_xy_1=[temp1_xy;ones(1,M_des2)];
                T_temp1_xy=solution*temp1_xy_1;
                temp2_xy=loc2(:,[1,2]);
                temp2_xy=temp2_xy';
                temp2_xy_2=[temp2_xy;ones(1,M_des2)];
                diff_xy=T_temp1_xy-temp2_xy_2;
                T_error=sqrt(sum(diff_xy.^2,1));

                JD=(1+T_error).*distance(k,:);

                [vals,index]=sort(JD);
                if(vals(1)/vals(2)<Td)
                    match(k,1)=index(1);
                    Dk=vals(1);
                    match(k,2)=Dk;
                else
                    match(k,1)=0;
                    match(k,2)=0;
                end

                if(i==1&&j==1)
                    match(k,:)=match(k,:);
                elseif(prev_match(k,1)==0)
                    match(k,:)=match(k,:);
                elseif(prev_match(k,1)~=0)
                    if(match(k,1)~=0&&match(k,2)<prev_match(k,2))
                        match(k,:)=match(k,:);
                    else
                        match(k,:)=prev_match(k,:);
                    end
                end
            end
        end
    end


    temp_match=match(:,1)';
    num=sum(temp_match>0);
    [~,point1,point2]=find(temp_match);
    loc11=loc1(point1,[1,2,3,4,5,6]);
    loc22=loc2(point2,[1,2,3,4,5,6]);
    loc11=[loc11,point2'];
    loc22=[loc22,point2'];


    [loc11,loc22,]=logic_filter(loc11,loc22,delat_x,delat_y,unit_x,unit_y,S,RMO);
    uni1=[loc11(:,[1,2]),loc22(:,[1,2])];
    [~,i,~]=unique(uni1,'rows','first');
    loc11=loc11(sort(i)',:);loc22=loc22(sort(i)',:);

    [end_solution,rmse,loc11,loc22]=FSC(loc11,loc22,change_form,3);
    point_1=loc11;
    point_2=loc22;

end