clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
sub_fir6 = {'S01','','S02','','S03','','S04','','S05','','S06'};
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
time=[0,100];
color = {[51, 34, 136],[17, 119, 51],[68, 170, 153],[136, 204, 238],[221, 204, 119],[204, 102, 119],[170, 68, 153],[136, 34 85]};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
%color = {'-r','-y','-g','-c','-b','-k'};
x = (1 : 2500) - 1300;
figure
for kk = 1 : 8
    
    sub_d_small = [];
    sub_d_large = [];
    
    for i = 1 : 2
        if i == 1 & kk > 6 
            continue
        end
        for j = 1 : 8

            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{kk}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
%            sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');

            ms_path= fullfile(MATpath,sprintf('%s.mat', direction_name{j}) );
            load(ms_path);

            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            
            timepoint = countfigure_new(MS_TEMP,tab,samplingRateData);
            
            small_o = tab(:,11) == 0.5;
            large_o = tab(:,11) == 8;
            
            timepoint_small = timepoint(small_o,:);
            timepoint_large = timepoint(large_o,:);
            
            s_mean_small = nanmean(timepoint_small,1)*60;
            s_mean_large = nanmean(timepoint_large,1)*60;
            
            sub_d_small = [sub_d_small;s_mean_small(1000:2194)];
            sub_d_large = [sub_d_large;s_mean_large(1000:2194)];
            
        end
    end
    subplot(2,1,1)
    hold on
    a = shadedErrorBar(x(1000:2194), mean(sub_d_small,1),std(sub_d_small,0,1)/sqrt(size(sub_d_small,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
    hold on
    b = plot(x(1000:2194), mean(sub_d_small,1),'-','LineWidth',1);
    b.Color = color{kk}/255;
    xline(0,'-','Stimuli on')
    hold on 
    xline(500,'-','Stimuli off')
    ylim([0,6])
    xlim([-300,894])

    ylabel('Microsaccade Rate (hz)')
    xlabel('Time (ms)')
    title('Smallest tile')
    subplot(2,1,2)
    hold on
    a1 = shadedErrorBar(x(1000:2194), mean(sub_d_large,1),std(sub_d_large,0,1)/sqrt(size(sub_d_large,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
    hold on
    b1 = plot(x(1000:2194), mean(sub_d_large,1),'-','LineWidth',1);
    b1.Color = color{kk}/255;
    xline(0,'-','Stimuli on')
    hold on 
    xline(500,'-','Stimuli off')
    ylim([0,6])
    xlim([-300,894])

    ylabel('Microsaccade Rate (hz)')
    xlabel('Time (ms)')
    title('Largest tile')
end
figpath='C:\Users\86186\Desktop\fig\new\time_series_SvL.png';
%title_file = [figpath,'/',subject{kk},'- within stimuli.png'];
saveas(gca,figpath)