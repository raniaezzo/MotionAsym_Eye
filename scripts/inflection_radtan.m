%% inflection for radial and tangential conditions
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'UR','LL','UL','LR','VU','VL','HR','HL'};
padding_time=[0,100];
c={'Switched','Original'};

for kk = 1 : 6          
    for i = 2 : 2
        for j = 1 : 8
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
            if length(tab)~=800
                a = find(tab(:,9)==0);
                dd=setdiff(1:length(tab),a);
                tab = tab(dd',:);
            end
            event_name = [subject{kk},'_',c{1},'_',direction{j},'_events.mat'];
            load(event_name)
            
            if j ==1 | j == 2 
                
                rad = EVENTS( tab(:,9) == 1 | tab(:,9) == 2,:);
                tan = EVENTS( tab(:,9) == 3 | tab(:,9) == 4,:);

            elseif j ==3 | j == 4
                tan = EVENTS( tab(:,9) == 1 | tab(:,9) == 2,:);
                rad = EVENTS( tab(:,9) == 3 | tab(:,9) == 4,:);
                                
            elseif j ==5 | j == 6
                tan = EVENTS( tab(:,9) == 7 | tab(:,9) == 8,:);
                rad = EVENTS( tab(:,9) == 5 | tab(:,9) == 6,:);                
               
            elseif j ==7 | j == 8
                rad = EVENTS( tab(:,9) == 7 | tab(:,9) == 8,:);
                tan = EVENTS( tab(:,9) == 5 | tab(:,9) == 6,:);    
            end 
            
            rate_rad = nansum(rad,1)./sum(~isnan(rad),1); 
            rate_tan = nansum(tan,1)./sum(~isnan(tan),1); 
            rate_name_rad = [subject{kk},'_',c{i},'_',direction{j},'_rate_rad.mat'];
            rate_name_tan = [subject{kk},'_',c{i},'_',direction{j},'_rate_tan.mat'];
            save(rate_name_rad,'rate_rad')
            save(rate_name_tan,'rate_tan')
            
        end
    end
end