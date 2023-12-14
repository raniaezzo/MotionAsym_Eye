clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
sub_fir6 = {'S01','','S02','','S03','','S04','','S05','','S06'};
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
time=[0,100];
color = {[51, 34, 136],[17, 119, 51],[68, 170, 153],[136, 204, 238],[221, 204, 119],[204, 102, 119],[170, 68, 153],[136, 34 85]};

%color = {'-r','-y','-g','-c','-b','-k'};
x = 1 : 2500;

ssd = [];
for kk = 1 : 6
    
    diff = [];
    
    for i = 1 : 2
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

            ms_path= fullfile(MATpath,sprintf('%s.mat', direction{j}) );
            load(ms_path);

            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            dd = tab(:,8) - tab(:,2);
            diff = [diff;dd];
        end
    end
    
    diff = diff(diff<3000);
    diff = diff(diff>0);
    ssd = [ssd ; mean(diff)];
end
            
            
            
            
            
            
            
            