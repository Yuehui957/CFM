function[solution,rmse,cor1,cor2,diagnostic]=match_ablation(des1,loc1,des2,loc2,change_form,gate_mode)
% Diagnostic variant of the supplied matching code; detector/descriptors unchanged.
% gate_mode: original | bypass | model_box (default).
% Requires the user's original ransac, FSC, LSM, XY_position_offset, logic_filter.
% Output columns retain the original match.m final order: [x y col5 col6 col3 col4 targetID].
    if nargin<6, gate_mode='model_box'; end
    gate_mode=char(gate_mode);
    assert(any(strcmp(gate_mode,{'original','bypass','model_box'})), 'Unknown gate mode.');
    assert(any(strcmp(change_form,{'Similarity','Affine'})), 'Use Similarity or Affine.');
    assert(size(loc1,2)>=6 && size(loc2,2)>=6, 'loc needs at least six columns.');
    assert(size(des1,1)==size(loc1,1) && size(des2,1)==size(loc2,1), 'Row count mismatch.');
    assert(size(des1,2)==size(des2,2) && size(des2,1)>=2, 'Invalid descriptor arrays.');
    assert(isreal(des1) && isreal(des2) && all(isfinite(des1(:))) && all(isfinite(des2(:))), 'Invalid descriptors.');
    diagnostic=struct('gate_mode',gate_mode);
    nmin=2+strcmp(change_form,'Affine');

    

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


    assert(size(cor1,1)>nmin, 'Too few candidates for the legacy sampler.');
    assert(size(unique(cor1(1:end-1,1:2),'rows'),1)>=nmin && size(unique(cor2(1:end-1,1:2),'rows'),1)>=nmin, 'Degenerate candidates for legacy sampler.');
    diagnostic.initial_candidates=[cor1(:,1:2),cor2(:,1:2)];
    [solution_old,~,cor1,cor2]=ransac(cor1,cor2,change_form,6);

    diagnostic.initial_inliers=[cor1(:,1:2),cor2(:,1:2)];
    diagnostic.initial_solution=solution_old;
    assert(all(isfinite(solution_old(:))) && size(cor1,1)>=nmin, 'Initial model failed.');
    [solution,cor1,cor2,rmse,stage]=cfm_ablation_refine(solution_old,...
    loc1,loc2,cor1,cor2,des1,des2,distRatio,change_form,cos_dis,gate_mode);
    diagnostic.stage=stage;

end




%%


function[end_solution,point_1,point_2,rmse,stage]=...
    cfm_ablation_refine(solution,loc1,loc2,cor1,cor2,des1,des2,Td,change_form,distance,gate_mode)

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


    stage.before_gate=[loc11(:,1:2),loc22(:,1:2)];
    stage.histogram=[S(1),RMO(1),delat_x,delat_y];
    switch gate_mode
        case 'original'
            [loc11,loc22]=logic_filter(loc11,loc22,delat_x,delat_y,unit_x,unit_y,S,RMO);
        case 'bypass'
            % Only remove the histogram gate; final FSC is still applied.
        case 'model_box'
            % Same rectangular pixel tolerances, using the existing full model.
            keep=cfm_model_box_mask(loc11(:,1:2),loc22(:,1:2),solution,unit_x,unit_y);
            loc11=loc11(keep,:); loc22=loc22(keep,:);
    end
    stage.after_gate=[loc11(:,1:2),loc22(:,1:2)];
    uni1=[loc11(:,[1,2]),loc22(:,[1,2])];
    [~,i,~]=unique(uni1,'rows','first');
    loc11=loc11(sort(i)',:);loc22=loc22(sort(i)',:);

    nmin=2+strcmp(change_form,'Affine');
    if size(loc11,1)<=nmin
        end_solution=nan(3); rmse=NaN;
        point_1=zeros(0,7); point_2=zeros(0,7);
        stage.final=zeros(0,4); stage.status='insufficient_candidates'; return;
    end
    [end_solution,rmse,loc11,loc22]=FSC(loc11,loc22,change_form,3);
    stage.final=[loc11(:,1:2),loc22(:,1:2)];
    stage.status='completed';
    point_1=loc11;
    point_2=loc22;

end


%%

function results=compare_cfm_matching(des1,loc1,des2,loc2,change_form,save_path)
% Reuse the SAME extracted features in all three arms. No detector reruns.
% Saves all stages as [source_x source_y target_x target_y].
    if nargin<6, save_path='cfm_matching_comparison.mat'; end
    modes={'original','bypass','model_box'};
    results=struct();
    fprintf('Mode          Candidates Initial Rematch Gate Final FitRMSE\n');
    for k=1:numel(modes)
        mode=modes{k}; rng(0,'twister');
        [H,r,c1,c2,d]=match_ablation(des1,loc1,des2,loc2,change_form,mode);
        results.(mode)=struct('H',H,'fit_rmse',r,'cor1',c1,'cor2',c2,'diagnostic',d);
        fprintf('%-13s %8d %7d %7d %4d %5d %.4f\n',mode,...
            size(d.initial_candidates,1),size(d.initial_inliers,1),...
            size(d.stage.before_gate,1),size(d.stage.after_gate,1),size(c1,1),r);
        if k>1
            assert(isequaln(d.initial_candidates,results.original.diagnostic.initial_candidates), 'Initial candidates differ.');
            assert(isequaln(d.initial_solution,results.original.diagnostic.initial_solution), 'Initial models differ.');
            assert(isequaln(d.stage.before_gate,results.original.diagnostic.stage.before_gate), 'Pre-gate candidates differ.');
        end
    end
    save(save_path,'results');
    figure('Name','CFM matching ablation: source-coordinate coverage');
    for k=1:numel(modes)
        d=results.(modes{k}).diagnostic;
        sets={d.stage.before_gate,d.stage.after_gate,d.stage.final};
        labels={'Before gate','After gate','After FSC'};
        for j=1:3
            subplot(3,3,(k-1)*3+j);
            plot(loc1(:,1),loc1(:,2),'.','Color',[.85 .85 .85]); hold on;
            v=sets{j}; plot(v(:,1),v(:,2),'r.');
            axis equal; set(gca,'YDir','reverse');
            title(sprintf('%s | %s | n=%d',modes{k},labels{j},size(v,1)), 'Interpreter','none');
            xlabel('Source x'); ylabel('Source y');
        end
    end
end


%%

function keep=cfm_model_box_mask(p,q,H,tx,ty)
% H maps source [x;y;1] to target coordinates; tolerances are target pixels.
    z=H*[p';ones(1,size(p,1))];
    predicted=bsxfun(@rdivide,z(1:2,:),z(3,:))';
    e=q-predicted;
    keep=all(isfinite(e),2) & abs(e(:,1))<tx & abs(e(:,2))<ty;
end




%%
function m=cfm_known_transform_error(cor1,cor2,H_est,H_true,check_xy,tolerance)
% Independent truth only. H_true maps image1 to image2 in pixel coordinates.
% check_xy: fixed Nx2 evaluation positions within the common valid region.
% Use the SAME positions/truth tolerance for all arms, not fitted inlier points.
    if nargin<6, tolerance=3; end
    assert(isequal(size(H_true),[3 3]) && size(check_xy,2)==2 && ~isempty(check_xy), 'Invalid truth or evaluation positions.');
    m=struct('n_matches',size(cor1,1),'n_correct',0,'precision',NaN,...
             'truth_match_rmse',NaN,'transform_rmse',NaN,'truth_errors',[],'correct_mask',false(size(cor1,1),1));
    if ~isempty(cor1)
        e=sqrt(sum((warp_points(H_true,cor1(:,1:2))-cor2(:,1:2)).^2,2));
        m.truth_errors=e; m.correct_mask=isfinite(e)&e<tolerance;
        m.n_correct=sum(m.correct_mask); m.precision=mean(m.correct_mask);
        m.truth_match_rmse=sqrt(mean(e.^2));
    end
    if isequal(size(H_est),[3 3]) && all(isfinite(H_est(:)))
        e=warp_points(H_est,check_xy)-warp_points(H_true,check_xy);
        m.transform_rmse=sqrt(mean(sum(e.^2,2)));
    end
end
function q=warp_points(H,p)
    z=H*[p';ones(1,size(p,1))];
    q=bsxfun(@rdivide,z(1:2,:),z(3,:))';
end
