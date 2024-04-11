%% this function to add outside num clone into tab file
function add_outside(tab_path,path_dat, samplingRateData)

    outside_path = strrep(path_dat, 'Dat_stim', 'outside');

    samplesPerMilisec = samplingRateData/1000;

    [path,filename, ~] = fileparts(tab_path);
    load(tab_path);
    load(outside_path);
    len=size(tab,1);
    bl=zeros(len,1);

    stimStart = samplesPerMilisec*1300;
    stimEnd = samplesPerMilisec*1800;
    for i = 1 : len
        order=tab(i,2)+stimStart <=outside(:,2) & tab(i,2)+stimEnd >=outside(:,2);
        c1=outside(:,1);
        num=sum(c1(order));
        bl(i)=num;
    end
    tab=[tab,bl];
%     save(sprintf('%s/%s_outside.mat',path,filename),'tab');
%     new_tab_outside=sprintf('%s/%s_outside.mat',path,filename);
    save(sprintf('%s/%s.mat',path,filename),'tab');
end