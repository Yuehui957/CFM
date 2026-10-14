function result = CFM_Fig4(gradient_space, K_abs, keypoint, out_dir, curvature_mode)
% CFM_Fig4 Export four spatially corresponding descriptor intermediate maps.
% Example (inside GLOH_descriptors_improve, before the for loop):
%   idx = 1;
%   CFM_Fig4(gradient,K_abs,key_point_array(idx,:), ...
%       fullfile(pwd,'CFM_Fig4_output'),'as_provided');
%
% keypoint = [x, y, scale, layer, main_angle]; angle is not needed here.
% curvature_mode:
%   'as_provided' (default): reproduce the supplied descriptor's indexing.
%       For 3-D K_abs, its two-subscript ROI accesses the first plane.
%   'matched_layer': use K_abs(:,:,layer). This is a deliberate change from
%       the supplied 3-D-array behavior; update the matcher consistently.
% No new gradients or curvatures are calculated: existing arrays are used.
% Images are display renderings; maps_raw.mat preserves numerical values.

if nargin < 4 || isempty(out_dir), out_dir=fullfile(pwd,'CFM_Fig4_output'); end
if nargin < 5 || isempty(curvature_mode), curvature_mode='as_provided'; end
validateattributes(gradient_space,{'numeric'},{'real','nonempty','finite','nonnegative'});
validateattributes(K_abs,{'numeric'},{'real','nonempty','finite'});
validateattributes(keypoint,{'numeric'},{'real','vector','finite','nonempty'});
if numel(keypoint)<4, error('keypoint needs [x,y,scale,layer].'); end
if ndims(gradient_space)>3 || ndims(K_abs)>3, error('Only 2-D/3-D arrays are supported.'); end
x=double(keypoint(1)); y=double(keypoint(2)); scale=double(keypoint(3)); layer=double(keypoint(4));
[H,W,ng]=size(gradient_space);
if H~=size(K_abs,1) || W~=size(K_abs,2), error('Gradient and curvature spatial sizes differ.'); end
if x~=round(x)||y~=round(y)||layer~=round(layer)||x<1||x>W||y<1||y>H||layer<1||scale<=0
    error('x,y,layer must be valid positive integer indices; scale must be positive.');
end
if ng==1
    g=double(gradient_space); % Explicit 2-D selected-layer input is supported.
else
    if layer>ng, error('Gradient layer is out of bounds.'); end
    g=double(gradient_space(:,:,layer));
end
nk=size(K_abs,3);
switch lower(char(curvature_mode))
    case 'as_provided'
        curvature_layer=1;
    case 'matched_layer'
        if nk==1, curvature_layer=1;
        else, curvature_layer=layer; end
        if curvature_layer>nk, error('Curvature layer is out of bounds.'); end
    otherwise
        error('curvature_mode must be as_provided or matched_layer.');
end
k=abs(double(K_abs(:,:,curvature_layer)));
radius=round(min(12*scale,min(H,W)/3));
x1=max(1,x-radius); x2=min(W,x+radius);
y1=max(1,y-radius); y2=min(H,y+radius);
A=g(y1:y2,x1:x2); sub_K=k(y1:y2,x1:x2);
maxK=max(sub_K(:));
if maxK>0, Knorm=sub_K/maxK; else, Knorm=zeros(size(sub_K)); end

% Exactly reproduce the supplied ROI-local 9x9 neighborhood calculation.
% Exclude the center; at ROI boundaries count valid neighbors only.
B=zeros(size(A)); half_win=4;
[nr,nc]=size(A);
for row=1:nr
    for col=1:nc
        total=0; count=0;
        for di=-half_win:half_win
            for dj=-half_win:half_win
                if di==0 && dj==0, continue; end
                rr=row+di; cc=col+dj;
                if rr>=1 && rr<=nr && cc>=1 && cc<=nc
                    total=total+abs(A(rr,cc)-A(row,col));
                    count=count+1;
                end
            end
        end
        if count>0, B(row,col)=total/count; else, B(row,col)=A(row,col); end
    end
end
alpha=0.2;
C=alpha+(1-alpha)*Knorm; % Actual range: [0.2,1], not [0,1].
D=B.*C;
if ~exist(out_dir,'dir'), mkdir(out_dir); end
out_dir=char(out_dir);
% A has its own scale. B and D share one scale, avoiding misleading rescaling.
maxA=max(A(:)); if maxA==0, maxA=1; end
maxBD=max([B(:);D(:)]); if maxBD==0, maxBD=1; end
limits=[0,maxA;0,maxBD;0,1;0,maxBD];
names={'a_gradient_magnitude','b_local_gradient_difference', ...
       'c_curvature_weight','d_weighted_magnitude'};
labels={'(a) Gradient magnitude','(b) Local gradient difference', ...
        '(c) Curvature weight','(d) Curvature-weighted magnitude'};
maps={A,B,C,D};
for n=1:4
    v=maps{n}/limits(n,2); v=min(max(v,0),1);
    % 16-bit grayscale preserves display detail; MAT remains authoritative.
    imwrite(uint16(round(v*65535)),fullfile(out_dir,[names{n},'.png']));
end
f=figure('Color','w','Name','CFM Figure 4','Position',[80,80,1050,850]);
for n=1:4
    ax=subplot(2,2,n,'Parent',f);
    imagesc(ax,maps{n},limits(n,:)); axis(ax,'image'); axis(ax,'off');
    colormap(ax,gray(256)); colorbar(ax);
    title(ax,labels{n},'FontName','Times New Roman','FontSize',13, ...
        'FontWeight','normal','Color','k','Interpreter','none');
end
annotation(f,'textbox',[0.04,0.002,0.94,0.04], ...
    'String','Same ROI in all panels. (b) and (d) share a display scale; (c) uses [0,1].', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',10,'Color','k');
set(f,'PaperPositionMode','auto','InvertHardcopy','off');
print(f,fullfile(out_dir,'Fig4_combined.png'),'-dpng','-r600');
print(f,fullfile(out_dir,'Fig4_combined.pdf'),'-dpdf','-bestfit','-painters');
savefig(f,fullfile(out_dir,'Fig4_editable.fig'));
result=struct('gradient_magnitude',A,'gradient_difference',B, ...
    'curvature_weight',C,'weighted_magnitude',D,'curvature_abs',sub_K, ...
    'keypoint',keypoint,'roi_xyxy',[x1,y1,x2,y2], ...
    'radius',radius,'gradient_layer',layer,'curvature_layer',curvature_layer, ...
    'curvature_mode',char(curvature_mode),'alpha',alpha,'half_window',half_win, ...
    'display_limits',limits);
save(fullfile(out_dir,'maps_raw.mat'),'result');
fid=fopen(fullfile(out_dir,'map_metadata.txt'),'w');
if fid>=0
    fprintf(fid,'Mode: %s\nGradient layer: %d\nCurvature array layer: %d\n', ...
        char(curvature_mode),layer,curvature_layer);
    fprintf(fid,'ROI [x1 y1 x2 y2]: %d %d %d %d\n',x1,y1,x2,y2);
    fprintf(fid,'alpha=0.2; 9x9 ROI-local window; center excluded.\n');
    fprintf(fid,'Each grayscale PNG is scaled by its display limits below.\n');
    for n=1:4
        fprintf(fid,'%s: raw range [%.16g %.16g], display range [%.16g %.16g]\n', ...
            names{n},min(maps{n}(:)),max(maps{n}(:)),limits(n,1),limits(n,2));
    end
    fclose(fid);
end
fprintf('Figure 4 exported to: %s\n',out_dir);
fprintf('Gradient layer=%d; curvature array layer=%d; mode=%s\n', ...
    layer,curvature_layer,char(curvature_mode));
end
