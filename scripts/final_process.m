clc; clear all;
homedir = '/Users/rania/Downloads/MS_Project';
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
condition_t = {'Switch','Original'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
time=[0,100];
vel = []; vel2 = []; vel3=[];
vel_sti = [];
am = []; am2 = []; am3=[];
am_sti = [];
for kk = 1 : 8 %8   % subject       
    padding_time=[0,100];
    for i = 1 : 2 %2 %2 % orig vs switch
        if i == 1 & kk > 6 
            continue
        end
        for j = 1 : 8 %4 %8 %dir

            main_folder = fullfile(homedir,'Data_DI_wEYE', subject{kk}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
%            sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');

            ms_path= fullfile(MATpath,sprintf('%s.mat', direction{j}) );
            load(ms_path);

            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            MS_TEMP = filterAM(MS_TEMP,[0.05,1]);
            save(sprintf('%s/%s.mat',MATpath, direction_name{j}),'MS_TEMP');
            
            AM=sqrt(MS_TEMP(:,6).^2+MS_TEMP(:,7).^2); %add total amplitude
            MS_TEMP=[MS_TEMP,AM];
            MS_TEMP=calculate_degree(MS_TEMP);
            MS_TEMP(MS_TEMP(:,12)<0,12)=MS_TEMP(MS_TEMP(:,12)<0,12)+360;
            MS_TEMP(MS_TEMP(:,13)<0,13)=MS_TEMP(MS_TEMP(:,13)<0,13)+360;
            
            if samplingRateData == 1000
                vel = [vel;MS_TEMP(:,3)];
                am = [am ; MS_TEMP(:,11)];
                disp('1000')
            elseif samplingRateData == 2000
                vel2 = [vel2;MS_TEMP(:,3)];
                am2 = [am2 ; MS_TEMP(:,11)];
                disp('2000')
            elseif samplingRateData == 500
                vel3 = [vel3;MS_TEMP(:,3)];
                am3 = [am3 ; MS_TEMP(:,11)];
                disp('500')
            end
            

            %find MS trial
            whole_trial = 1 : size(MS_TEMP,1);
            lower=(1300-time(1))*samplingRateData/1000;
            higher=(1800+time(2))*samplingRateData/1000;
            ms_start=MS_TEMP(:,8);
            ms_end=MS_TEMP(:,9);
%            ms_trial=MS_TEMP(:,10);
            a= lower<ms_start & higher > ms_start;
            trial_w1=whole_trial(a);
            b= lower<ms_end & higher > ms_end;
            trial_w2=whole_trial(b);
            trial_w=unique([trial_w1,trial_w2]); 
            
            vel_sti = [vel_sti;MS_TEMP(trial_w,3)];
            am_sti = [am_sti;MS_TEMP(trial_w,11)];
 
            
%             figure
%             scatter(am_sti,vel_sti)
%             title([subject{kk},'-',condition{i},'-',direction{j},'- within stimuli'])
%             ylabel('amplitude')
%             xlabel('max velocity')
% 
%             figpath='C:\Users\86186\Desktop\fig';
%             title_file = [figpath,'/',subject{kk},'- within stimuli.png'];
            %saveas(gca,title_file)

%             figure
%             scatter(am,vel)
%             ylabel('amplitude')
%             xlabel('max velocity')
%             title([subject{kk},'/',condition{i},'/',direction{j}])
%             figpath='C:\Users\86186\Desktop\fig';
%             title_file_all = [figpath,'/',subject{kk},'- all.png'];
            %saveas(gca,title_file_all)

            
            
        end
%         figure
%         scatter(am,vel)
%         ylabel('amplitude')
%         xlabel('max velocity')
%         xlim([0,1])
%         ylim([0,80])
%         title([subject{kk},'/',condition_t{i}])
%         figpath='C:\Users\86186\Desktop\fig';
%         title_file_all = [figpath,'/',subject{kk},'-',condition_t{i},'- all.png'];
%         saveas(gca,title_file_all)
    end
%     figure
%     scatter(am_sti,vel_sti)
%     title([subject{kk},'- within stimuli'])
%     ylabel('amplitude')
%     xlabel('max velocity')
% 
%     figpath='C:\Users\86186\Desktop\fig';
%     title_file = [figpath,'/',subject{kk},'- within stimuli.png'];
%     saveas(gca,title_file)
%             
%     figure
%     scatter(am,vel)
%     ylabel('amplitude')
%     xlabel('max velocity')
%     title(subject{kk})
%     figpath='C:\Users\86186\Desktop\fig';
%     title_file_all = [figpath,'/',subject{kk},'- all.png'];
%     saveas(gca,title_file_all)
end

%% 
figure
s = scatter(am,vel, 'k', 'filled', 'MarkerFaceAlpha', .02);
hold on
s = scatter(am2,vel2/2, 'k', 'filled', 'MarkerFaceAlpha', .02);
hold on
s = scatter(am3,vel3*2, 'k', 'filled', 'MarkerFaceAlpha', .02);
ylabel('Peak Velocity (deg/s)')
xlabel('Amplitude (deg)')
xlim([0,1])
ylim([0,80])
%title('Main Sequence')
s.SizeData = 35;
ylim([0 60])
xlim([0.05 1])
ax = gca;
ax.YTick = [0 15 30 45 60];
ax.XTick = [0.2 0.4 0.6 0.8 1];
ax.FontSize = 18;
f1 = gcf;
f1.Position = [864 479 560 202];
%figpath='C:\Users\86186\Desktop\fig';
%title_file_all = 'C:\Users\86186\Desktop\fig\All_AvV.png';
%saveas(gca,title_file_all)

%%
figure
mdl = fitlm(am_sti,vel_sti);
plot(mdl);

legend('off')
ylabel('Peak Velocity (deg/s)')
xlabel('Amplitude (deg)')
xticks([0.05,0.2:0.1:1])
xlim([0.05,1])
ylim([0,80])
title('Main Sequence')
%figpath='C:\Users\86186\Desktop\fig';
%title_file_all = 'C:\Users\86186\Desktop\fig\new\All_AvV_sti.png';
%saveas(gca,title_file_all)


%mdl_sti = fitlm(am_sti,vel_sti);