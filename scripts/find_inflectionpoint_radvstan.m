%% find inflection point for radial and tang
clear all; clc;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';
c={'Switched','Original'};
Trial_w = [];
Time_point = zeros(6,2,8,3);
Value = zeros(6,2,8,3);
for ii = 1:6
%     if ii == 5
%         continue
%     end
    for jj = 2 : 2
        for kk = 1 : 8
            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{ii}, ...
                    'RawData', condition{jj}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            rate_name = [subject{ii},'_',c{jj},'_',direction{kk},'_rate_tan.mat'];
            edf_name = dir(sprintf('*%s*.edf', direction{jj})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
%             samplingRateData=findSamplingRate(msg_filepath);
%             x=1:1000/samplingRateData:2500;
            load(rate_name)
%             if ii ==5
            rate1=rate_tan(1:length(rate_tan)*0.8);
%             end
            titlename = [subject{ii},' ',c{jj},' ',direction{kk}];
            time=[500,2000];
            pa_pass=0.001;
            temporalthresh=70;
            onset=1;
            stepsize=1;
%             if ii == 5                 
            [t1,t2,t3,y_fil] = infl(rate1,time,pa_pass,temporalthresh,onset,stepsize);
%             else
%                 [t1,t2,t3,y_fil] = infl(rate_rad,time,pa_pass,temporalthresh,onset,stepsize);
%             end
            if isempty(t3)
                t3 = 2000*length(rate_tan)/2500;
                trial = [ii,jj,kk];
                Trial_w = [Trial_w;trial];
            end
%             x=1:length(rate_rad);
%             x=x./(length(rate_rad)/2500);
%             if ii == 5                 
%                 y_fil = [y_fil,rate_rad(length(rate)*0.8+1:length(rate))];
%             end
            v = rate_tan([t1,t2,t3]);
            t=[t1,t2,t3]*2500/length(rate_tan);
            Time_point(ii,jj,kk,:) = t;
            Value(ii,jj,kk,:) = v;
            save('Point_Value_tan.mat','Value')
            save('Point_time_tan.mat','Time_point')
            save('No_inflection_tan.mat','Trial_w')
%% for draw figure
%             figure
%             subplot(2,1,1)
%             plot(x,rate)
%             subplot(2,1,2)
%             plot(x,y_fil)
%             hold on 
%             xline(t1/(length(rate)/2500))
%             hold on 
%             xline(t2/(length(rate)/2500))
%             hold on 
%             xline(t3/(length(rate)/2500))
%             title(titlename)
%             figpath='C:\Users\86186\Desktop\fig';
%             name = [figpath,'/',subject{ii},'_',c{jj},'_',direction{kk},'_p.png'];
%             saveas(gca,name)            
                
                
                
                
                
                
                
        end
    end
end
save('Point_Value_tan.mat','Value')
save('Point_time_tan.mat','Time_point')
save('No_inflection_tan.mat','Trial_w')

