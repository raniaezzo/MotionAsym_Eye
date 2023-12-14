%% calculate the suppression and inflection point
clear all; clc;
homedir = '/Users/rania/Downloads/MS_Project';
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
%edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';
c={'Switched','Original'};
%% Switched vs. Original Condition
load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\Point_Value.mat');
load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\Point_time.mat');
rad = load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\S06\RawData\Full_distance_radialtangential\Block1\eyedata\Point_time_rad.mat');
tan = load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\S06\RawData\Full_distance_radialtangential\Block1\eyedata\Point_time_tan.mat');
y1=ones(8,1);
y2=ones(8,1)*2;
y3=ones(8,1)*3;
y4=ones(8,1)*4;
y5=ones(8,1)*5;
y6=ones(8,1)*6;
y_mean=[1 2 3 4 5 6];
rad_time = rad.Time_point;
tan_time = tan.Time_point;


%% plot
for i = 1 : 6 
    sub_pol_car = squeeze(Time_point(i,2,:,3));
    sub_pol_obl = squeeze(Time_point(i,1,:,3));
    sub_car_car = squeeze(Time_point(i,:,1:4,3));
    sub_car_car = reshape(sub_car_car,[8 1]);
    sub_car_obl = squeeze(Time_point(i,:,5:8,3));
    sub_car_obl = reshape(sub_car_obl,[8 1]);
    sub_rad = squeeze(rad_time(i,2,:,3));
    sub_tan = squeeze(tan_time(i,2,:,3));
    mean_point = squeeze(mean([sub_pol_car,sub_pol_obl,sub_car_car,sub_car_obl,sub_rad,sub_tan],1));
    figure
    scatter(mean_point,y_mean,'r','filled')
    hold on
    scatter(sub_pol_car,y1,'b')
    hold on 
    scatter(sub_pol_obl,y2,'b')
    hold on 
    scatter(sub_car_car,y3,'b')
    hold on 
    scatter(sub_car_obl,y4,'b')
    hold on 
    scatter(sub_rad,y5,'b')
    hold on 
    scatter(sub_tan,y6,'b')
    hold on 
    xline(1300)
    hold on 
    xline(1800)
    legend('mean')
    title(subject{i})
    xlabel('time(ms)')
    xlim([1200 2100])
    ylim([0 7])
    yticks(y_mean)
    yticklabels({'Polar Cardinal','Polar Oblique','Cartesion Cardinal','Cartesion Oblique','Radial','Tangential'})
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/',subject{i},'_inflection.png'];
    saveas(gca,name)  
end