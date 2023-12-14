% plot correlation
clear all; clc;
x = 1 : 8 ; 
subject = {'S01','S02','S03','S04','S05','S06'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
path = 'F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE\MS_rate_and_acc_correlation.mat';
load(path)
c = squeeze(mean(correlation,1));
sem = zeros(2,8);
for j = 1 : 8 
    sem(1,j) = std(squeeze(correlation(:,1,j)))/sqrt(6);
    sem(2,j) = std(squeeze(correlation(:,1,j)))/sqrt(6);
end



figure
for i = 1 : 6 
    subplot(5,2,i)
    scatter(x,squeeze(correlation(i,1,:)),'b')
    hold on 
    scatter(x,squeeze(correlation(i,2,:)),'r')
    hold on
    xticks(x)
    xticklabels(direction)
    title(subject{i})
end
% legend({'Switched Condition','Original Condition'},'Location','east');
% legend('boxoff')
subplot(5,2,[7,8,9,10])
scatter(x,c(1,:)','b')
hold on 
scatter(x,c(2,:)','r')
hold on
errorbar(x,c(1,:)',sem(1,:)','bo')
hold on 
errorbar(x,c(2,:)',sem(2,:)','ro')
xlim([0,9])
xticks(x)
xticklabels(direction)
title('ALL Subject')
legend({'Switched Condition','Original Condition'},'Location','east');
legend('boxoff')
%legend({'Switched Condition','Original Condition'},'Location','east');
figpath='C:\Users\86186\Desktop\fig';
name = [figpath,'/','subject_co.png'];
saveas(gca,name)
%lgd.Layout.Tile = 9;