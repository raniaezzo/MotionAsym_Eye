%% plot the raw data figure
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
con = {'Location Oblique','Location Cardinal'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_s3 = {'VL','HL','HR','LL','LR','UL','UR'};
le_direction = {'VU','VL','HL','HR','LL','LR','UL','UR','mean-car','mean-obl'};
le_direction_s3 = {'VL','HL','HR','LL','LR','UL','UR','mean-car','mean-obl'};

loc_direction ={'Cardinal Direction','Oblique Direction'};

tilt_direction = {'0.5°','1°','2°','4°','8°'};

color={[215,25,28],[253,174,97],[255,255,191],[171,217,233],[44,123,182]};	

num = [];
diff=[];
x=(1:2201) - 1300;
figpath='C:\Users\86186\Desktop\fig';



pupil_data1 = []; %0.5
pupil_data2 = []; %1
pupil_data3 = []; %2
pupil_data4 = []; %4
pupil_data5 = []; %8
%extract pupil size data
for kk = 7:8
   pupil_data11 = []; %0.5
   pupil_data12 = []; %1
   pupil_data13 = []; %2
   pupil_data14 = []; %4
   pupil_data15 = []; %8
%     if kk == 3
%          continue
%     end
    
    pupil_data = [];

    for i = 1 :2
    
        if kk > 6 & i == 1 
            continue
        end
                
        for j = 1 : 8
            if i == 1 & j ==1 & kk ==3 
                continue
            end
            pupil_data = [];
            main_folder = fullfile('F:\pupildata\Data_DI_wEYE\Data_DI_wEYE', subject{kk}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
    %        sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            cd('MATs')
            cdt=dir([direction{j},'*_Dat_all_new_22.mat']);
            load(cdt.name)
            
            for kkk = 1 : 800
                nondata = 0 ;
                trial_st = tab(kkk,2) ;
                aa=Dat_all_new(trial_st + 500<=Dat_all_new(:,1) & Dat_all_new(:,1)<=trial_st+2200,:);
                if isempty(aa) ==1 
                    nondata = nondata +1 ;
                    continue
                end
                num = [num,size(aa,1)];
%                 diff = [diff,max(aa(:,1))-min(aa(:,1))];
                [~,I]= sort(aa(:,1));
                timepoint = aa(I,4)';
                
                
                if tab(kkk,11) == 0.5                    
                    if samplingRateData == 500
                        timepoint = upsample_raw(timepoint);
                        if size(timepoint,2) == 1699
                            timepoint = [timepoint(1),timepoint,timepoint(length(timepoint))];                
                        end
                    end
                    if samplingRateData == 2000
                        timepoint = downsample(timepoint,2);
                    end
                    if size(timepoint,2) == 1701  
                        pupil_data11 = [pupil_data11;timepoint];
                    end 
                elseif tab(kkk,11) == 1
                    if samplingRateData == 500
                        timepoint = upsample_raw(timepoint);
                        if size(timepoint,2) == 1699
                            timepoint = [timepoint(1),timepoint,timepoint(length(timepoint))];                
                        end
                    end
                    if samplingRateData == 2000
                        timepoint = downsample(timepoint,2);
                    end
                    if size(timepoint,2) == 1701  
                        pupil_data12 = [pupil_data12;timepoint];
                    end
                elseif tab(kkk,11) == 2
                    if samplingRateData == 500
                        timepoint = upsample_raw(timepoint);
                        if size(timepoint,2) == 1699
                            timepoint = [timepoint(1),timepoint,timepoint(length(timepoint))];                
                        end
                    end
                    if samplingRateData == 2000
                        timepoint = downsample(timepoint,2);
                    end
                    if size(timepoint,2) == 1701  
                        pupil_data13 = [pupil_data13;timepoint];
                    end
                elseif tab(kkk,11) == 4
                    if samplingRateData == 500
                        timepoint = upsample_raw(timepoint);
                        if size(timepoint,2) == 1699
                            timepoint = [timepoint(1),timepoint,timepoint(length(timepoint))];                
                        end
                    end
                    if samplingRateData == 2000
                        timepoint = downsample(timepoint,2);
                    end
                    if size(timepoint,2) == 1701  
                        pupil_data14 = [pupil_data14;timepoint];
                    end
                else 
                    if samplingRateData == 500
                        timepoint = upsample_raw(timepoint);
                        if size(timepoint,2) == 1699
                            timepoint = [timepoint(1),timepoint,timepoint(length(timepoint))];                
                        end
                    end
                    if samplingRateData == 2000
                        timepoint = downsample(timepoint,2);
                    end
                    if size(timepoint,2) == 1701  
                        pupil_data15 = [pupil_data15;timepoint];
                    end
                    
                    
                end
             end

            

        end

    end
    pupil_data11=pupil_data11(pupil_data11(:,1701)~=0,:);
    pupil_data12=pupil_data12(pupil_data12(:,1701)~=0,:);
    pupil_data13=pupil_data13(pupil_data13(:,1701)~=0,:);
    pupil_data14=pupil_data14(pupil_data14(:,1701)~=0,:);
    pupil_data15=pupil_data15(pupil_data15(:,1701)~=0,:);
    baseline1  = mean(pupil_data11(:,600:800),2);
    baseline2  = mean(pupil_data12(:,600:800),2);
    baseline3  = mean(pupil_data13(:,600:800),2);
    baseline4  = mean(pupil_data14(:,600:800),2);
    baseline5  = mean(pupil_data15(:,600:800),2);
    
    pupil_data11 = (pupil_data11-baseline1)./baseline1*100; 
    pupil_data12 = (pupil_data12-baseline2)./baseline2*100; 
    pupil_data13 = (pupil_data13-baseline3)./baseline3*100; 
    pupil_data14 = (pupil_data14-baseline4)./baseline4*100; 
    pupil_data15 = (pupil_data15-baseline5)./baseline5*100; 

    pupil_data1 = [pupil_data1;mean(pupil_data11,1)];
    pupil_data2 = [pupil_data2;mean(pupil_data12,1)];
    pupil_data3 = [pupil_data3;mean(pupil_data13,1)];
    pupil_data4 = [pupil_data4;mean(pupil_data14,1)];
    pupil_data5 = [pupil_data5;mean(pupil_data15,1)];



end


figure
hold on 
shadedErrorBar(x(1001:2201), mean(pupil_data1(:,501:1701),1),std(pupil_data1(:,501:1701),0,1)/sqrt(size(pupil_data1(:,501:1701),1)),'lineprops',{'-','Color',color{1}/255},'transparent',1,'patchSaturation',0.05);
hold on 
shadedErrorBar(x(1001:2201), mean(pupil_data2(:,501:1701),1),std(pupil_data2(:,501:1701),0,1)/sqrt(size(pupil_data2(:,501:1701),1)),'lineprops',{'-','Color',color{2}/255},'transparent',1,'patchSaturation',0.05);
hold on 
shadedErrorBar(x(1001:2201), mean(pupil_data3(:,501:1701),1),std(pupil_data3(:,501:1701),0,1)/sqrt(size(pupil_data3(:,501:1701),1)),'lineprops',{'-','Color',color{3}/255},'transparent',1,'patchSaturation',0.05);
hold on 
shadedErrorBar(x(1001:2201), mean(pupil_data4(:,501:1701),1),std(pupil_data4(:,501:1701),0,1)/sqrt(size(pupil_data4(:,501:1701),1)),'lineprops',{'-','Color',color{4}/255},'transparent',1,'patchSaturation',0.05);
hold on 
shadedErrorBar(x(1001:2201), mean(pupil_data5(:,501:1701),1),std(pupil_data5(:,501:1701),0,1)/sqrt(size(pupil_data5(:,501:1701),1)),'lineprops',{'-','Color',color{5}/255},'transparent',1,'patchSaturation',0.05);
hold on 
ab1=plot(x(1001:2201), mean(pupil_data1(:,501:1701),1),'-b','LineWidth',2);
hold on
ab2=plot(x(1001:2201), mean(pupil_data2(:,501:1701),1),'-r','LineWidth',2);
hold on 
ab3=plot(x(1001:2201), mean(pupil_data3(:,501:1701),1),'-g','LineWidth',2);
hold on
ab4=plot(x(1001:2201), mean(pupil_data4(:,501:1701),1),'-k','LineWidth',2);
hold on
ab5=plot(x(1001:2201), mean(pupil_data5(:,501:1701),1),'-y','LineWidth',2);
ab1.Color = color{1}/255;
ab2.Color = color{2}/255;
ab3.Color = color{3}/255;
ab4.Color = color{4}/255;
ab5.Color = color{5}/255;


a1=xline(-100,'-','Baseline','alpha',0.065);
a=xline(-200:0,'-y','alpha',0.065);
a2 = xline(0,'-' ,'Stimuli on');
a3 = xline(500,'-','Stimuli off');
legend(tilt_direction{:},'Location','northeastoutside')
xlabel('Time(ms)')
ylabel('Pupil area (% change from baseline)') 
title_name = 'All Subject';
%         figpath='C:\Users\86186\Desktop\fig';
name = 'C:\Users\86186\Desktop\fig\all_data_tilt.png';
title(title_name)
saveas(gca,name)

