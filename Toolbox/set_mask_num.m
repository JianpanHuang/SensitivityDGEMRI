function num = set_mask_num(varargin)
% Jianpan Huang 20190630 
% Set the number of masks that will ben drawn

    dlg_titl = ['Set mask number'];
    prom = {['Mask number is:']};
    num_line = [1 50];
    defa_answ = {'1'};
    options.Resize ='on';
    answ = inputdlg(prom, dlg_titl, num_line, defa_answ, options);
    num = str2double(answ{1});
end