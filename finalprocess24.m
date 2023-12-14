% find main reaction time for per subject
clc; clear all;

homedir = '/Users/rania/Downloads/MS_Project';
tab_total = [] ;
MS_total = [] ; 
sample=[];
sub={'S01','S02','S03','S04','S05','S06','S07','S08'};
subject = 'S01'; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};

figpath='C:\Users\86186\Desktop\fig\new\new_pre';
color = {[127, 191, 123]	,[175, 141, 195]}	;

% color={[215,25,28],[253,174,97],[255,255,191],[171,217,233],[44,123,182]};	
react = zeros([1,8]);

for k = 1 : 8
    rec_sub = [];
    for i = 1 : 2
        if i ==1 & k > 6
            continue
        end
        for j = 1 : 8

            main_folder = fullfile(homedir,'Data_DI_wEYE', sub{k}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
            sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
    %         cd(fullfile(main_folder, 'eyedata','MATs'));
    %         ms_name = dir(sprintf('%s.mat', direction{j})).name;
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction_name{j}) );
            load(ms_path);
    %         MS_TEMP(:,10)=MS_TEMP(:,10)+ ((i-1)*8 + (j-1))*800;   
    %         MS_total=[MS_total;MS_TEMP];
    %         MATpath = fullfile(main_folder, 'eyedata','MATs');
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            rec_sub = [rec_sub;tab(:,14)];
        end
    end
    react(k) = mean(rec_sub(rec_sub < 1.5))*1000;
end


name='reaction.mat';
save(fullfile(homedir,name),'react')    
            