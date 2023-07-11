clear all;clc; close all;
addpath(genpath('toolbox'))
%% Parameters
data_path = 'Data/C57mouse_50%glc1';
exp_num = 24;
mask_name = 'csf_small';
method_num = 2; % 1 for T2, 2 for VDMP, 3 for T1roh
method_name = {'CPMG','onVDMP','onSL'};
base_num = 8;
time_resol = 90 ; % Time resolution, in s
disca_num = 2; % Discard the first N images because of non-steady-state
delet_num = 0;
base_time = base_num*time_resol/60;
csf_thresh = 0.1752; % A ratio between 0~1
if_denoi = 0; % 0 without denoising, 1 with denoising 

%% Load image data
data_dir = [data_path, filesep, int2str(exp_num), filesep, 'Result_2dseq.mat'];
load(data_dir);
img = Result.image;
% img(:,:,1:6) are the images for 
% T2 Par, T2 CSF, T2 VDMP, T2 VDMP, T1roh Par, T1roh CSF, respectively (then repeated)
img(:,:,1:disca_num*6)=[];
img_paren = img(:,:,method_num:6:end);
img_csf = img(:,:,3+method_num:6:end);

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
[u, s] = mlsvd(img_csf, [med_ind(1), med_ind(2), nel_ind(3)]);
if if_denoi == 1
    img_csf = lmlragen(u, s);
end
[xs, ys, ts] = size(img_csf);

%% Draw ROI
[mask, mask_num] = draw_mask_csf(data_path, img_csf(:,:,1), [mask_name,'.mat'], 'gray', csf_thresh);

%% Calculate DGE
img_base = mean(img_csf(:,:,1:base_num),3);
for mm = 1:ts
    img_dge(:,:,mm)=(img_base - img_csf(:,:,mm))./img_base;
end

%% DGE point-wise calculation
dge_sig = zeros(ts,mask_num);
for m = 1:ts
    img_temp = img_dge(:,:,m);
	for n = 1:mask_num
        roi_temp = mask(:,:,n);
        dge_sig(m,n) = mean2(img_temp(roi_temp==1)); 
	end
end
time_min = (time_resol*(1:ts))'/60; % Time in minutes
dge_sig = dge_sig*100; % Transfer to percentage

%% Fitting and display
t0 = base_time + time_resol/60;
%     Smax     miu_out  miu_in
p0 = [1        1        1  ];
lb = [1        0        0];
ub = [100      100      100];
time_min_cut = time_min(1:end-delet_num);
s = @(p,t) p(1)*(exp(-p(2)*(t-t0))./(1+exp(-p(3)*(t-t0))));
% plusfun = @(x) max(x,0);
% s = @(p,t)  p(1)-p(1)*exp(-plusfun(t-p(2))/p(3));
for n = 1:mask_num
    s_raw = dge_sig(1:end-delet_num, n);
    [p_fit(n,:), rn(n,1)] = lsqcurvefit(s, p0, time_min_cut, s_raw, lb, ub);
    s_fit(:,n) = s(p_fit(n,:), time_min);
    r2(n,1) = 1 - rn(n,1)/sum(s_raw.^2);
end
varia = var(dge_sig-s_fit);
snrc = 10*log10(p_fit(:,1)^2./varia);

% Display
yscal = [-5, 25];
xscal = [(1-base_num)*time_resol/60 max(time_min)-base_time];
figure
subplot(1,2,1)
imagesc(img_csf(:,:,1)); axis image; colormap(gray); hold on;
contour(mask,1,'m-','LineWidth',2); title('ROI'); set(gca,'Fontsize',16);
subplot(1,2,2)
for n = 1:mask_num
    plot(time_min-base_time,dge_sig(:,n),'bo',time_min-base_time,s_fit(:,n),'b-','LineWidth',2)
    hold on
end 
xlabel('Time (min)','FontName', 'Arial','FontSize',16);
ylabel('Normalized Intensity (%)','FontName', 'Arial','FontSize',16);
title(['DGE cureve (',mask_name,')'])
axis([xscal yscal]);
set(gcf,'Position',[300 300 1300 500]);
set(gca,'Fontsize',16)
% Show the fit parameters on figure
txt_str = ['Smax = ',num2str(p_fit(1)),', uin = ',num2str(p_fit(2)),', SNR = ',num2str(snrc)];
txt_ui = uicontrol('Style', 'text',...
                  'Fontname', 'Arial',...
                  'FontSize', 16,...
                  'Position', [748 420 400 30],...
                  'String', txt_str);
hold off
% Save DGE results
saveas(gcf,[data_path,filesep,num2str(exp_num),filesep,'CSF_DGE_',cell2mat(method_name(method_num)),'.jpg']);
if if_denoi == 1
    xlswrite( [data_path,filesep,int2str(exp_num),filesep,'CSF_DGE_',cell2mat(method_name(method_num)),'_denoi.xls'],[time_min-base_time,dge_sig,time_min-base_time,s_fit]);
else
    xlswrite( [data_path,filesep,int2str(exp_num),filesep,'CSF_DGE_',cell2mat(method_name(method_num)),'_raw.xls'],[time_min-base_time,dge_sig,time_min-base_time,s_fit]);
end
% Save fit parameters
smax = p_fit(:,1);
uin = p_fit(:,2);
if if_denoi ==1
    save_path = [data_path,filesep,num2str(exp_num),filesep,'CSF_Smax_uin_R2_SNR_',cell2mat(method_name(method_num)),'_denoi.txt'];
else
    save_path = [data_path,filesep,num2str(exp_num),filesep,'CSF_Smax_uin_R2_SNR_',cell2mat(method_name(method_num)),'_raw.txt'];
end
save_mat = [smax uin r2 snrc];
save_txt(save_path,save_mat);