%% this function add blink into tab file
function add_blink(tab_path,blink_path,samplingRateData)

    samplesPerMilisec = samplingRateData/1000;

    [path,filename, ~] = fileparts(tab_path);
    load(tab_path)
    load(blink_path)
    tab_len=size(tab,1);
    c1=zeros(tab_len,1);

    stimStart = samplesPerMilisec*1300;
    stimEnd = samplesPerMilisec*1800;

    for i = 1 : tab_len
        order=tab(i,2)+stimStart <=blink(:,1) & tab(i,2)+stimEnd >=blink(:,1);
        c1(i)=sum(order);
    end
    tab=[tab,c1];
    %save(sprintf('%s/%s_blink.mat',path,filename),'tab');
    %path=sprintf('%s/%s_blink.mat',path,filename);
    save(sprintf('%s/%s.mat',path,filename),'tab');
end