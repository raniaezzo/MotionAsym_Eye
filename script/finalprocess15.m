clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
location = {'UR','LL','UL','LR','VU','VL','HR','HL'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
oblique_cum = zeros(5,8);
cardinal_cum = zeros(5,8);
tilt_angle = [0.5,1,2,4,8];
time = [0,100];
for kk = 7 : 8 
    sz = [320 9];
    varNames = {'Subject', 'tilt angle', 'motion direction', 'location', 'amplitude (degrees)', 'peak velocity (degrees)', 'duration (ms)', 'total # trials with at least 1 MS', 'total # MS'};
    varTypes = {'string','double','string','string','double','double','double','double','double'};
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

    T.Subject(:) = subject{kk};
        
    padding_time=[0,100];
    for i = 1 : 2
        if kk > 6 & i == 1 
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
            AM = sqrt(MS_TEMP(:,6).^2 + MS_TEMP(:,7).^2);
            MS_TEMP = [MS_TEMP,AM];
            start_or = (i-1) * 160 + (j -1 ) * 20+1;
            end_or  = (i-1) * 160 + j * 20;
            
            
            T(start_or:end_or,3) = direction(j);
            loc_vec = tab(:,9);
            tilt_vec = tab(:,11);
            for loc = 1 : 4 
                if tab(1,9) > 4 
                    loc_s = loc + 4;
                else
                    loc_s = loc ;
                end 
                for til = 1:5                    
                    T((start_or-1+(loc - 1 )*5 + til),4)=location(loc_s);
                    T((start_or-1+(loc - 1 )*5 + til),2)={tilt_angle(til)};
                    index_loc = find(loc_vec==loc_s);
                    index_tilt = find(tilt_vec==tilt_angle(til));
                    final_index = intersect(index_loc,index_tilt);
                    ms_index = [];
                    for tt = 1 : length(final_index) 
                        ms_index = [ms_index;find(MS_TEMP(:,10)==final_index(tt))];
                    end
                    MS_index = MS_TEMP(ms_index,:);
                    whole_trial = 1 : size(MS_index,1);
                    lower=(1300-time(1))*samplingRateData/1000;
                    higher=(1800+time(2))*samplingRateData/1000;
                    ms_start=MS_index(:,8);
                    ms_end=MS_index(:,9);
                    a= lower<ms_start & higher > ms_start;
                    trial_w1=whole_trial(a);
                    b= lower<ms_end & higher > ms_end;
                    trial_w2=whole_trial(b);
                    trial_w=unique([trial_w1,trial_w2]); 
                    if isempty(trial_w)
                        continue
                    end
                    MS_within = MS_index(trial_w,:);
                    
                    T((start_or-1+(loc - 1 )*5 + til),5)={mean(MS_within(:,11))};%am
                    T((start_or-1+(loc - 1 )*5 + til),6)={mean(MS_within(:,3))};% peak velo
                    duration = MS_within(:,2)-MS_within(:,1);
                    T((start_or-1+(loc - 1 )*5 + til),7)={mean(duration)};%duration
                    T((start_or-1+(loc - 1 )*5 + til),9)={size(MS_within,1)}; % num
                    T((start_or-1+(loc - 1 )*5 + til),8)={size(unique(MS_within(:,10)),1)};
%                    T = T(161:320,:);

                end
            end                             
        end
    end
    figpath = 'C:\Users\86186\Desktop\fig\new\table';
    name = [figpath,'/',subject{kk},'-t.mat'];
    save(name,'T')
end
