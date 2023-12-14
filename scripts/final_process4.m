%% plot the raw data figure
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
subject_le = {'S01','S02','S03','S04','S05','S06'};
subject_le2 = {'S01','','S02','','S03','','S04','','S05','','S06','','S07','','S08',''};
con = {'Location Oblique','Location Cardinal'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_s3 = {'VL','HL','HR','LL','LR','UL','UR'};
le_direction = {'VU','VL','HL','HR','LL','LR','UL','UR','mean-car','mean-obl'};
le_direction_s3 = {'VL','HL','HR','LL','LR','UL','UR','mean-car','mean-obl'};

loc_direction ={'Cardinal Direction','Oblique Direction'};

%color = {'-r','-y','-g','-c','-b','-k'};

color = {[51, 34, 136],[17, 119, 51],[68, 170, 153],[136, 204, 238],[221, 204, 119],[204, 102, 119],[170, 68, 153],[136, 34 85]};
num = [];
diff=[];
% color = {'-r','-b','-g','-y','--r','--b','--g','--y'};
%color = {'-r','-b'};
x=(1:2201) - 1300;
figpath='C:\Users\86186\Desktop\fig';
%extract pupil size data
figure
for kk = 1:8
        if kk == 5 
            continue
        end
        if kk<7
            meanfile1 = [figpath,'/',subject{kk},'-CL-CD-22.mat'];
            meanfile2 = [figpath,'/',subject{kk},'-CL-OD-22.mat'];
            load(meanfile1)
            load(meanfile2)           
            baseline  = mean(pupil_data_1(:,600:800),2);
            pupil_data_1 = (pupil_data_1-baseline)./baseline*100;
            baseline  = mean(pupil_data_2(:,600:800),2);
            pupil_data_2 = (pupil_data_2-baseline)./baseline*100;           
            meanfile3 = [figpath,'/',subject{kk},'-OL-CD-22.mat'];
            meanfile4 = [figpath,'/',subject{kk},'-OL-OD-22.mat'];
            load(meanfile3)
            load(meanfile4)           
            baseline  = mean(pupil_data_3(:,600:800),2);
            pupil_data_3 = (pupil_data_3-baseline)./baseline*100;
            baseline  = mean(pupil_data_4(:,600:800),2);
            pupil_data_4 = (pupil_data_4-baseline)./baseline*100; 
            
            pupil_data = [pupil_data_1;pupil_data_2;pupil_data_3;pupil_data_4];
        else
            meanfile1 = [figpath,'/',subject{kk},'-CL-CD-22.mat'];
            meanfile2 = [figpath,'/',subject{kk},'-CL-OD-22.mat'];
            load(meanfile1)
            load(meanfile2)           
            baseline  = mean(pupil_data_1(:,600:800),2);
            pupil_data_1 = (pupil_data_1-baseline)./baseline*100;
            baseline  = mean(pupil_data_2(:,600:800),2);
            pupil_data_2 = (pupil_data_2-baseline)./baseline*100;           
            
            pupil_data = [pupil_data_1;pupil_data_2];
        end     
            
        hold on
        ab1 = plot(x(1001:2201), mean(pupil_data(:,501:1701),1),'-','LineWidth',3);
        hold on
        ab2 = shadedErrorBar(x(1001:2201), mean(pupil_data(:,501:1701),1),std(pupil_data(:,501:1701),0,1)/sqrt(size(pupil_data(:,501:1701),1)),'lineprops',{'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.05);
        hold on 
        ab1.Color = color{kk}/255;
        %ab2.Color = color{kk}/255;
end

a1=xline(-100,'-','Baseline','alpha',0.065);
a=xline(-200:0,'-y','alpha',0.065);
a2 = xline(0,'-','Stimuli on');
a3 = xline(500,'-','Stimuli off');
%legend(subject_le2,'Location','northeastoutside')
xlabel('Time(ms)')
ylabel('Pupil area (% change from baseline)')
xlim([-300,900])
%title_name = [subject{kk},'- ',con{i},];

name = 'C:\Users\86186\Desktop\fig\new\subject.png';
%title('All subject')
saveas(gca,name) 
