%% for matlab figure to png
path = 'C:\Users\86186\Desktop\fig';
filelist = dir(path);

for i = 1:1:length(filelist)
    suf = strsplit(filelist(i).name, '.');
    if length(suf) < 2 % 不是以后缀名结尾的文件
        contine
    else
        if strcmp(suf{2}, 'fig') == 1
            fig = openfig(strcat(path, '\', filelist(i).name));
            saveas(fig, strcat(path, '\', suf{1}, '.png'), 'png') %保存为png格式的图片
        end    
    end
end
