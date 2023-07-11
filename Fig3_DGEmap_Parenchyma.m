clear all;clc; close all;
addpath(genpath('toolbox'))

%% Parameters
data_path = 'Data/C57mouse_50%glc1';
exp_num = 24;
mask_name = 'brain';
method_num = 2; % 1 for T2, 2 for VDMP, 3 for T1roh
method_name = {'CPMG', 'onVDMP', 'onSL'};
base_num = 8;
time_resol = 90 ; % Time resolution, in s
disca_num = 2; % Discard the first N images because of non-steady-state
delet_num = 0;
base_time = base_num*time_resol/60;
save_path = 'Fig/Fig3';

%% Load image data
data_dir = [data_path, filesep, int2str(exp_num), filesep, 'Result_2dseq.mat'];
load(data_dir);
img = Result.image;
% img(:,:,1:6) are the images for 
% T2 Par, T2 CSF, T2 VDMP, T2 VDMP, T1roh Par, T1roh CSF, respectively (then repeated)
img(:,:,1:disca_num*6)=[];
img_paren = img(:,:,method_num:6:end);

%% Denoising
[~, ~, sv] = mlsvd(img_paren); % Singular value
svn{1} = sv{1}/max(sv{1}); % Normalized singular value
svn{2} = sv{2}/max(sv{2}); 
svn{3} = sv{3}/max(sv{3}); 
% Determine the truncation indexes using Malinowskis, Nelson and Median criteria.
[mal_ind(1,1), nel_ind(1,1), med_ind(1,1)] = trunc_determ(svn{1});
[mal_ind(1,2), nel_ind(1,2), med_ind(1,2)] = trunc_determ(svn{2});
[mal_ind(1,3), nel_ind(1,3), med_ind(1,3)] = trunc_determ(svn{3});
% Denoise
[u, s] = mlsvd(img_paren, [med_ind(1), med_ind(2), nel_ind(3)]);
img_paren = lmlragen(u, s);
[xs, ys, ts] = size(img_paren);

%% Draw ROI
[mask, mask_num] = draw_mask(data_path, img_paren(:,:,1), [mask_name,'.mat'], 'gray');

%% Calculate DGE
img_base = mean(img_paren(:,:,1:base_num),3);
for mm = 1:ts
    img_dge(:,:,mm)=(img_base - img_paren(:,:,mm))./img_base;
end
img_dge = img_dge.*mask*100;

%% Display
sample_num = 4;
for n = 1:floor(ts/sample_num)
    img_dge_sample(:,:,n) = mean(img_dge(:,:,(n-1)*sample_num+1:n*sample_num),3);
end
% xscl = 19:54; 
xscl = 19:60; 
yscl = 20:75; 
img_dge_sample_brain = img_dge_sample(xscl,yscl,:);
% lastDgeImg = sum(dgeImgEnh(:,:,152-7:152),3)/smpNum*100;
fst_img = zeros(length(xscl),length(yscl));
img_dge_sample_brain = cat(3,fst_img,img_dge_sample_brain);
img_dge_sample_brain_2d = imshow3dimage(img_dge_sample_brain, 4);
h = figure('numbertitle','off','name','DGEImg','color','white');
imshow(img_dge_sample_brain_2d,[],'InitialMagnification','fit');setsaveas(h,600, 600);
mycolormap(1); colorbar; caxis([0 5]); colorbar('FontWeight','bold','FontName','Arial','FontSize',12); drawnow;
% set(gcf,'Position',[100 300 1200 800]);
export_fig([save_path, filesep, 'Parenchyma_DGEmap_every',num2str(sample_num),'_',cell2mat(method_name(method_num))], '-jpg', '-r200')
anatom_img = img_paren(:,:,1).*mask;
anatom_img_save = anatom_img(xscl,yscl);
imwrite(anatom_img_save/max(anatom_img_save(:)),[save_path,filesep,'Parenchyma_anatomy','_',cell2mat(method_name(method_num)),'.tiff'],'tiff');

%% Last DGE map for Graphical Abstract
h2 = figure('numbertitle','off','name','LastDGEImg','color','white');
imshow(img_dge_sample(xscl,yscl,end),[],'InitialMagnification','fit');setsaveas(h2,600, 600);
mycolormap(1); caxis([0,5]); drawnow;
export_fig([save_path, filesep, 'Last_Parenchyma_DGEmap_',cell2mat(method_name(method_num))], '-jpg', '-r200')