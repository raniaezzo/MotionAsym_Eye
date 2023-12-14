%% slope calculation and plot
clear all; clc;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';
c={'Switched','Original'};
load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\Point_Value.mat');
load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\Point_time.mat');
rad = load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\S06\RawData\Full_distance_radialtangential\Block1\eyedata\Point_time_rad.mat');
tan = load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\S06\RawData\Full_distance_radialtangential\Block1\eyedata\Point_time_tan.mat');
rad_v = load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\S06\RawData\Full_distance_radialtangential\Block1\eyedata\Point_Value_rad.mat');
tan_v = load('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\S06\RawData\Full_distance_radialtangential\Block1\eyedata\Point_Value_tan.mat');
%% calculate the slope
total_slope = zeros(6,2,8);
for i = 1: 6
    for j = 1 :2 
        for k = 1 : 8 
            total_slope(i,j,k)=(Value(i,j,k,2)-Value(i,j,k,1))/(Time_point(i,j,k,2)-Time_point(i,j,k,1));
        end
    end
end

rad_slope = zeros(6,8);
tan_slope = zeros(6,8);

for i =1 :6 
    for j = 1 : 8 
        rad_slope(i,j)=(rad_v.Value(i,2,j,2)-rad_v.Value(i,2,j,1))/(rad.Time_point(i,2,j,2)-rad.Time_point(i,2,j,1));
        tan_slope(i,j)=(tan_v.Value(i,2,j,2)-tan_v.Value(i,2,j,1))/(tan.Time_point(i,2,j,2)-tan.Time_point(i,2,j,1));
    end
end
        
%% plot
x1=ones(8,1);
x2=ones(8,1)*2;
x3=ones(8,1)*3;
x4=ones(8,1)*4;
x5=ones(8,1)*5;
x6=ones(8,1)*6;
x_mean=[1 2 3 4 5 6];



for i = 1 : 6 
    sub_pol_car = squeeze(total_slope(i,2,:))*1000;
    sub_pol_obl = squeeze(total_slope(i,1,:))*1000;
    sub_car_car = squeeze(total_slope(i,:,1:4))*1000;
    sub_car_car = reshape(sub_car_car,[8 1]);
    sub_car_obl = squeeze(total_slope(i,:,5:8))*1000;
    sub_car_obl = reshape(sub_car_obl,[8 1]);
    sub_rad = squeeze(rad_slope(i,:))'*1000;
    sub_tan = squeeze(tan_slope(i,:))'*1000;
    mean_point = squeeze(mean([sub_pol_car,sub_pol_obl,sub_car_car,sub_car_obl,sub_rad,sub_tan],1));
    figure
    scatter(x_mean,mean_point,'r','filled')
    hold on
    scatter(x1,sub_pol_car,'b')
    hold on 
    scatter(x2,sub_pol_obl,'b')
    hold on 
    scatter(x3,sub_car_car,'b')
    hold on 
    scatter(x4,sub_car_obl,'b')
    hold on 
    scatter(x5,sub_rad,'b')
    hold on 
    scatter(x6,sub_tan,'b')
    legend('mean')
    title(subject{i})
    xlabel('Condition')
    ylabel('Slope (MS_rate(%)/time(s))')
    xlim([0 7])
%     ylim([0 7])
    xticks(x_mean)
    xticklabels({'Polar Cardinal','Polar Oblique','Cartesion Cardinal','Cartesion Oblique','Radial','Tangential'})
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/',subject{i},'_slope.png'];
    saveas(gca,name)  
end
            