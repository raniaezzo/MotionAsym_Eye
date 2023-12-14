%this is withinsession script for draw figures
clear all; clc; 
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
window_size = 80;
step_size = 10;
padding_time=[0,100];
num_win = ceil((800-window_size)/step_size);
Acc = zeros(16,num_win+1);
Rate= zeros(16,num_win+1);
% correlation = zeros(6,2,8);
for kk = 1 : 6
    for i = 1 : 2
        for j = 1 : 8
            if kk == 3 & i== 1 & j ==1 
                continue
            end
            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{kk}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
            MATpath = fullfile(main_folder, 'eyedata','MATs');
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction{j}) );
            load(ms_path);
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            [acc, rate_trial,rate_sti]  = sw_withinsession(tab,MS_TEMP, window_size,samplingRateData,padding_time,step_size);
            x = 1 : length(acc);
            x= step_size * x + window_size/2;
            ms_rate = rate_sti/window_size;
            correlation(kk,i,j) = corr(acc,ms_rate);
            Acc((i-1)*8+j,:) = acc;
            Rate((i-1)*8+j,:) = ms_rate;
%             figure
%             plot(x,acc-mean(acc),'r')
%             hold on
%             plot(x,rate_sti/window_size-mean(rate_sti/window_size),'b')
           
%             ylim([-0.3,0.3])  
%             c={'original','switched'};
%             titlename = [subject{kk},c{i}];
%             title(titlename,'FontSize',12);
%             figpath='C:\Users\86186\Desktop\fig';
%             name = [figpath,'/',subject{kk},c{i},'_acc.png'];
            
            
            
        end

    
    end
    acc_m = mean(Acc,1);
    rate_m=mean(Rate,1);
    figure
    plot(x,acc_m,'r','LineWidth',2)
    hold on 
    plot(x,rate_m,'b','LineWidth',2)
    for ii = 1 : 16 
            a=plot(x,Acc(ii,:),'r');
            a.Color(4)=0.3;
            hold on 
            b=plot(x,Rate(ii,:),'b');
            b.Color(4)=0.3;
            hold on   
    end
    xlabel('Trials')
    ylabel('Percentage')
    legend({'accuracy','MS rate'},'Location','east')
    legend('boxoff')
    title(subject{kk})
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/',subject{kk},'_insess.png'];
    saveas(gca,name)


    
 
end

