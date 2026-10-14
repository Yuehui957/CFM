function[solution,rmse,cor1,cor2]=match(des1,loc1,des2,loc2,change_form)
    

   % Td=0.9;
    distRatio=0.9;
    des2t=des2';

    parfor i=1:size(des1,1)

        dotprods=des1(i,:)*des2t;
        similerDistance=acos(dotprods);
        
        [vals,indx]=sort(similerDistance);
        cos_dis(i,:)=similerDistance';

        if(vals(1)<distRatio*vals(2))
            match(i)=indx(1);
        else
            match(i)=0;
        end
    end


%%

    [~,point1,point2]=find(match);
    cor1=loc1(point1,[1,2,3,4,5,6]);
    cor2=loc2(point2,[1,2,3,4,5,6]);
    cor1=[cor1,point2'];cor2=[cor2,point2'];

    uni1=[cor1(:,[1,2]),cor2(:,[1,2])];
    [~,i,~]=unique(uni1,'rows','first');
    cor1=cor1(sort(i)',:);
    cor2=cor2(sort(i)',:);


    [solution_old,~,cor1,cor2]=ransac(cor1,cor2,change_form,6);

    [solution,cor1,cor2,rmse]=Horizontal_vertical_position_erro(solution_old,...
    loc1,loc2,cor1,cor2,des1,des2,distRatio,change_form,cos_dis);

end






