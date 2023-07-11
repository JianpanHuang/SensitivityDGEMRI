function [mask, num] = draw_mask(path, img, name, color_style)
% Draw mask(s) for a given image. 
% AUTHOR: Jianpan Huang, jianpanhuang@outlook.com
%
% INPUT:
%                path: Pathway where the image was saved
%                 img: Image used for mask drawing
%                name: Name of mask (string).
%         color_style: Color style of image displaying
%
% OUTPUT:
%                mask: Mask matrix
%                 num: Number of mask(s)

%%
if nargin < 4
    color_style = 'gray';
end

if ~exist([path, filesep, name],'file') ==1
    num = set_mask_num();
    [xs, ys] = size(img);
    mask = zeros(xs, ys, num);
    for mm = 1:num
        scrSize = get(0,'ScreenSize');
        h1 = figure;
        set(h1,'Position',[scrSize(3)*0.2 scrSize(4)*0.3 scrSize(3)*0.4 scrSize(4)*0.6]); 
        imagesc(img/max(img(:))),colormap(color_style);
        title(['Please draw the mask #',num2str(mm)]);
        mask(:,:,mm) = roipoly;
        % Confirm the mask
        close gcf;
        mask_temp = mask(:,:,mm);
        figure
        imagesc(abs(img/max(img(:)))); axis image; colormap(color_style);
        hold on
        contour(mask_temp, 1, 'm-', 'LineWidth',2);
        title('Mask contour')
        uicontrol('Position',[350 2 200 30], 'String','Continue', 'Callback','uiresume(gcbf)');
        uiwait(gcf);
        close gcf;
    end
    save([path, filesep, name], 'mask') 
else
    mask_load = load([path, filesep, name]); % Load mask
    mask_cell = struct2cell(mask_load);
    mask = cell2mat(mask_cell);
    mask=logical(mask);
    num = size(mask,3);
    % Display
    for mm = 1:num
        mask_temp = mask(:,:,mm);
        figure
        imagesc(abs(img/max(img(:)))); axis image; colormap(color_style)
        hold on
        contour(mask_temp, 1, 'm-', 'LineWidth', 2);
        title('Mask contour')
        uicontrol('Position',[350 2 200 30], 'String','Continue', 'Callback','uiresume(gcbf)');
        uiwait(gcf);
        close gcf;
    end
end