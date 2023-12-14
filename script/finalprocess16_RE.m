clc; clear all;

homedir = '/Users/rania/Desktop/Data_DI_wEYE_bogeng/Data_DI_wEYE/';

subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
sub_fir6 = {'S01','','S02','','S03','','S04','','S05','','S06'};
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
time=[0,100];
color = {[51, 34, 136],[17, 119, 51],[68, 170, 153],[136, 204, 238],[221, 204, 119],[204, 102, 119],[170, 68, 153],[136, 34 85]};
%direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
direction_name = {'VU','VL','HL','HR','LL','LR','UL','UR'};
%color = {'-r','-y','-g','-c','-b','-k'};
x = (1 : 2500) - 1300;
figure
for kk = 1 : 8
    
    sub_d_car = [];
    sub_d_obl = [];
    
    for i = 1 : 2
        if i == 1 & kk > 6 
            continue
        end
        for j = 1 : 8

            main_folder = fullfile(homedir, subject{kk}, ...
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
            s_mean = nanmean(timepoint,1)*60;
            if j < 5 
                sub_d_car = [sub_d_car;s_mean(1000:2194)];
            else
                sub_d_obl = [sub_d_obl;s_mean(1000:2194)];
            end
                   
%             sub_d = [sub_d;s_mean(1000:2194)];
            
        end
    end
    subplot(2,1,1)
    hold on
    a = shadedErrorBar(x(1000:2194), mean(sub_d_car,1),std(sub_d_car,0,1)/sqrt(size(sub_d_car,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
    hold on
    b = plot(x(1000:2194), mean(sub_d_car,1),'-','LineWidth',1);
    b.Color = color{kk}/255;
    xline(0,'-','Stimuli on')
    hold on 
    xline(500,'-','Stimuli off')
    ylim([0,6])
    xlim([-300,894])

    ylabel('Microsaccade Rate (hz)')
    xlabel('Time (ms)')
    title('Cardinal')
    subplot(2,1,2)
    hold on
    a1 = shadedErrorBar(x(1000:2194), mean(sub_d_obl,1),std(sub_d_obl,0,1)/sqrt(size(sub_d_obl,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
    hold on
    b1 = plot(x(1000:2194), mean(sub_d_obl,1),'-','LineWidth',1);
    b1.Color = color{kk}/255;
    xline(0,'-','Stimuli on')
    hold on 
    xline(500,'-','Stimuli off')
    ylim([0,6])
    xlim([-300,894])

    ylabel('Microsaccade Rate (hz)')
    xlabel('Time (ms)')
    title('Oblique')
end

% for hh = 1 : 6    
%     hold on
%     plot(x, sub_d(hh,:),color{hh},'LineWidth',1)
%     hold on
% end
% xline(0,'-','Stimuli on')
% hold on 
% xline(500,'-','Stimuli off')
% ylim([0,6])
% xlim([-300,894])
% 
% ylabel('Microsaccade Rate (hz)')
% xlabel('Time (ms)')



%figpath='C:\Users\86186\Desktop\fig\new\time_series_CvO.png';
%title_file = [figpath,'/',subject{kk},'- within stimuli.png'];
%saveas(gca,figpath)