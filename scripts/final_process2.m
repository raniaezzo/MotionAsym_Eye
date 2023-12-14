clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
time=[0,100];
for kk = 7 : 8       
    padding_time=[0,100];
    vel = [];
    vel_sti = [];
    am = [];
    am_sti = [];
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
            
            AM=sqrt(MS_TEMP(:,6).^2+MS_TEMP(:,7).^2); %add total amplitude
            MS_TEMP=[MS_TEMP,AM];
            MS_TEMP=calculate_degree(MS_TEMP);
            MS_TEMP(MS_TEMP(:,12)<0,12)=MS_TEMP(MS_TEMP(:,12)<0,12)+360;
            MS_TEMP(MS_TEMP(:,13)<0,13)=MS_TEMP(MS_TEMP(:,13)<0,13)+360;
            
            ms_add = [];
            for ss = MS_TEMP(:,10)
                ms_add = [ms_add;tab(ss,9:15)];
            end
            
            MS_new = [MS_TEMP,ms_add];
            figpath = 'C:\Users\86186\Desktop\fig\new\table';
            filename =  [direction{j},'-ms_new.mat'];
            save(filename,'MS_new')
            
            
            
            
        end
    end

end
