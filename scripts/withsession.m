function [acc, rate_trial,rate_sti] = withsession(tab,MS_TEMP,samplingRateData,time)
time_dur = 0 : 80 : 800;
acc=zeros(10,1);
rate_trial=zeros(10,1);
rate_sti=zeros(10,1);
lower=(1300-time(1))*samplingRateData/1000;
higher=(1800+time(2))*samplingRateData/1000;
ms_start=MS_TEMP(:,8);
ms_end=MS_TEMP(:,9);
ms_trial=MS_TEMP(:,10);
a= lower<ms_start & higher > ms_start;
trial_w1=ms_trial(a);
b= lower<ms_end & higher > ms_end;
trial_w2=ms_trial(b);
trial_w=unique([trial_w1;trial_w2]);
for i = 1 : 10
    acc(i) = sum(tab(time_dur(i)+1:time_dur(i+1),15))/80;
    MS=MS_TEMP(time_dur(i)<MS_TEMP(:,10)& MS_TEMP(:,10) <=time_dur(i+1),:);
    rate_trial(i) =  size(MS,1);
    rate_sti(i)=sum(time_dur(i)<trial_w & trial_w <=time_dur(i+1));
end
end