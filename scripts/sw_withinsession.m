%shift window within session
function  [acc, rate_trial,rate_sti]  = sw_withinsession(tab,MS_TEMP, window_size,samplingRateData,time,step_size)
num_win = ceil((800-window_size)/step_size);
start_time = 0 : step_size : step_size*num_win;
end_time = window_size : step_size : window_size + step_size*num_win;
if end_time(num_win+1) > 800
    end_time(num_win+1) = 800;
end
acc=zeros(num_win+1,1);
rate_trial=zeros(num_win+1,1);
rate_sti=zeros(num_win+1,1);
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

for i = 1 : num_win+1
    acc(i) = sum(tab(start_time(i)+1:end_time(i),15))/80;
    MS=MS_TEMP(start_time(i)<MS_TEMP(:,10)& MS_TEMP(:,10) <=end_time(i),:);
    rate_trial(i) =  size(MS,1);
    rate_sti(i)=sum(start_time(i)<trial_w & trial_w <=end_time(i));
end


end