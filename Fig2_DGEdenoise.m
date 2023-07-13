clear all;clc; close all;
addpath(genpath(pwd))

%% Set parameters
data_path = 'Data/C57mouse_50%glc1';
exp_num = 24;
roi_paren = 'brain';
roi_csf = 'csf';
disc_num = 2; % Discard the first N images because of non-steady-state
method_num = 1; % 1 for T2, 2 for VDMP, 3 for T1roh
method_name = {'T2','VDMP','T1roh'};
save_path = 'Fig/Fig2';

%% Load image data
data_dir = [data_path, filesep, num2str(exp_num), filesep, 'Result_2dseq.mat'];
load(data_dir);
img = Result.image;
% img(:,:,1:6) are the images for 
% T2 Par, T2 CSF, T2 VDMP, T2 VDMP, T1roh Par, T1roh CSF, respectively (then repeated)
img(:,:,1:disc_num*6)=[];
img_paren = img(:,:,1:6:end); % Parenchyma T2 meted
img_csf = img(:,:,5:6:end); % CSF VDMP
[xs, ys, ts] = size(img_paren); 

%% Draw masks
% Signal masks
[mask_paren, ~] = draw_mask(data_path, img_paren(:,:,1), [roi_paren,'.mat'], 'gray');
[mask_csf, ~] = draw_mask_csf(data_path, img_csf(:,:,1), [roi_csf,'.mat'], 'gray', 0.07);
% Noise mask
len = 15;
mask_noi = zeros(size(mask_paren));
mask_noi(2:len+1, 2:len+1) = 1;
mask_noi(2:len+1, end-len:end-1) = 1;
mask_noi(end-len:end-1, 2:len+1) = 1;
mask_noi(end-len:end-1, end-len:end-1) = 1;

%% Calculate the singular values and determine the truncation indexes
[~, ~, sv] = mlsvd(img_paren); % Singular value
svn{1} = sv{1}/max(sv{1}); % Normalized singular value
svn{2} = sv{2}/max(sv{2}); 
svn{3} = sv{3}/max(sv{3}); 
% Determine the truncation indexes using Malinowskis, Nelson and Median criteria.
[mal_ind(1,1), nel_ind(1,1), med_ind(1,1)] = trunc_determ(svn{1});
[mal_ind(1,2), nel_ind(1,2), med_ind(1,2)] = trunc_determ(svn{2});
[mal_ind(1,3), nel_ind(1,3), med_ind(1,3)] = trunc_determ(svn{3});

%% Display the results
scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[scrsz(3)*0.05, scrsz(4)*0.4, scrsz(3)*0.6, scrsz(4)*0.32]);
set(0,'defaultfigurecolor','w') 
% [ha,pos]=tight_subplot(Nh,Nw, gap, marg_h, marg_w)
ha = tight_subplot(1,3,[.09, .06],[.18, .08],[.07, .02]);
% Sigular value as a function of retained components
% X
axes(ha(1)), plot(1:xs,svn{1},'-b','LineWidth',2), hold on,
plot([mal_ind(1,1),mal_ind(1,1)],[0, 1],'-.m','LineWidth',2), hold on,
plot([nel_ind(1,1),nel_ind(1,1)],[0, 1],'-.g','LineWidth',2), hold on,
plot([med_ind(1,1),med_ind(1,1)],[0, 1],'-.r','LineWidth',2), hold off;
axis([0, xs, 0, 1]); set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',17, 'LineWidth', 1.5); 
ylabel('Normalized sigular value', 'FontName','Arial', 'FontWeight','bold', 'FontSize',18); 
xlabel('Component number', 'FontName', 'Arial', 'FontWeight','bold', 'FontSize',18);
title('x', 'FontName','Arial', 'FontWeight','bold', 'FontSize',18);
% Y
axes(ha(2)), plot(1:ys,svn{2},'-b','LineWidth',2), axis([0, ys, 0, 1]); hold on,
plot([mal_ind(1,2),mal_ind(1,2)],[0, 1],'-.m','LineWidth',2), hold on,
plot([nel_ind(1,2),nel_ind(1,2)],[0, 1],'-.g','LineWidth',2), hold on,
plot([med_ind(1,2),med_ind(1,2)],[0, 1],'-.r','LineWidth',2), hold off;
axis([0, ys, 0, 1]); set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',17, 'LineWidth', 1.5); 
xlabel('Component number', 'FontName', 'Arial', 'FontWeight','bold', 'FontSize',18);
title('y', 'FontName','Arial', 'FontWeight','bold', 'FontSize',18);
% T
axes(ha(3)), plot(1:ts,svn{3},'-b','LineWidth',2), axis([0, ts, 0, 1]); hold on,
plot([mal_ind(1,3),mal_ind(1,3)],[0, 1],'-.m','LineWidth',2), hold on,
plot([nel_ind(1,3),nel_ind(1,3)],[0, 1],'-.g','LineWidth',2), hold on,
plot([med_ind(1,3),med_ind(1,3)],[0, 1],'-.r','LineWidth',2), hold off;
legend('Sigular value','Malinowskis', 'Nelson', 'Median')
axis([0, ts, 0, 1]); set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',17, 'LineWidth', 1.5); 
xlabel('Component number', 'FontName', 'Arial', 'FontWeight','bold', 'FontSize',18);
title('n', 'FontName','Arial', 'FontWeight','bold', 'FontSize',18);

export_fig([save_path, filesep, 'MLSVD'], '-jpg', '-r200');


%% Image denoising comparison ?parenchyma?
[u, s] = mlsvd(img_paren, [med_ind(1), med_ind(2), nel_ind(3)]);
img_paren_denoi =lmlragen(u, s);
% Compare pre and post denoising images
img_paren_demo = img_paren(:,:,end);
img_paren_denoi_demo = img_paren_denoi(:,:,end);
img_paren_sig = mean2(img_paren_demo(mask_paren==1));
img_paren_snr = img_snr(img_paren_demo, mask_paren, mask_noi);
img_paren_denoi_sig = mean2(img_paren_denoi_demo(mask_paren==1));
img_paren_denoi_snr = img_snr(img_paren_denoi_demo, mask_paren, mask_noi);
% Display
figure,imshow(abs(img_paren_demo),[],'initialMagnification','fit');colormap(gray);colorbar,caxis([0,20])
hold on
contour(mask_paren, 1, 'r-', 'LineWidth', 2);
hold on
contour(mask_noi, 1, 'g-', 'LineWidth', 2);
set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',28); 
text(30, 75, ['Mean = ', num2str(roundn(img_paren_sig, -2)),'0'], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(24, 85, ['SNR = ',num2str(roundn(img_paren_snr, -2)), ' dB'], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(40, 5, 'Raw', 'Color','yellow', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 28)
export_fig([save_path, filesep, 'Paren_img_raw'], '-jpg', '-r200');
%
figure,imshow(abs(img_paren_denoi_demo),[],'initialMagnification','fit');colormap(gray);colorbar,caxis([0,20])
set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',28); 
text(30, 75, ['Mean = ', num2str(roundn(img_paren_denoi_sig, -2))], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(24, 85, ['SNR = ',num2str(roundn(img_paren_denoi_snr, -2)), ' dB'], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(31, 5, 'Denoised', 'Color','yellow', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 28)
export_fig([save_path, filesep, 'Paren_img_denoi'], '-jpg', '-r200')
%
figure,imshow(abs(img_paren_demo - img_paren_denoi_demo)./max(img_paren_demo(:))*100,[],'initialMagnification','fit'),mycolormap(1),colorbar,caxis([0,5])
set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',28); 
text(30, 5, 'Difference', 'Color','yellow', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 28)
export_fig([save_path, filesep, 'Paren_img_diff'], '-jpg', '-r200')

%% Image denoising comparison ?CSF?
[u, s] = mlsvd(img_csf, [med_ind(1), med_ind(2), nel_ind(3)]);
img_csf_denoi =lmlragen(u, s);
% Compare pre and post denoising images
img_csf_demo = img_csf(:,:,end);
img_csf_denoi_demo = img_csf_denoi(:,:,end);
img_csf_sig = mean2(img_csf_demo(mask_csf==1));
img_csf_snr = img_snr(img_csf_demo, mask_csf, mask_noi);
img_csf_denoi_sig = mean2(img_csf_denoi_demo(mask_csf==1));
img_csf_denoi_snr = img_snr(img_csf_denoi_demo, mask_csf, mask_noi);
% Display
figure,imshow(abs(img_csf_demo),[],'initialMagnification','fit'); colormap(gray); colorbar; caxis([0,10])
hold on
contour(mask_csf, 1, 'r-', 'LineWidth', 2);
hold on
contour(mask_noi, 1, 'g-', 'LineWidth', 2);
set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',28); 
text(30, 75, ['Mean = ', num2str(roundn(img_csf_sig, -2))], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(24, 85, ['SNR = ',num2str(roundn(img_csf_snr, -2)), ' dB'], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(40, 5, 'Raw', 'Color','yellow', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 28)
export_fig([save_path, filesep, 'CSF_img_raw'], '-jpg', '-r200');
%
figure,imshow(abs(img_csf_denoi_demo),[],'initialMagnification','fit');colormap(gray); colorbar; caxis([0,10])
set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',28); 
text(30, 75, ['Mean = ', num2str(roundn(img_csf_denoi_sig, -2))], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(24, 85, ['SNR = ',num2str(roundn(img_csf_denoi_snr, -2)), ' dB'], 'Color','white', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 24)
text(31, 5, 'Denoised', 'Color','yellow', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 28)
export_fig([save_path, filesep, 'CSF_img_denoi'], '-jpg', '-r200')
%
figure,imshow(abs(img_csf_demo - img_csf_denoi_demo)./max(img_csf_demo(:))*100,[],'initialMagnification','fit'),mycolormap(1),colorbar,caxis([0,5])
set(gca, 'FontName','Arial', 'FontWeight','bold', 'FontSize',28); 
text(30, 5, 'Difference', 'Color','yellow', 'FontName','Arial', 'FontWeight','bold', 'FontSize', 28)
export_fig([save_path, filesep, 'CSF_img_Diff'], '-jpg', '-r200')