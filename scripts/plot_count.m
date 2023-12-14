% plot figure
function plot_count(tab,MS_TEMP,samplerate,time)
lower=(1300-time)*samplerate/1000;
higher=(1800+time)*samplerate/1000;
ms_start=MS_TEMP(:,8);
ms_end=MS_TEMP(:,9);
ms_trial=MS_TEMP(:,10);
a= lower<ms_start & higher > ms_start;
trial_w1=ms_trial(a);
b= lower<ms_end & higher > ms_end;
trial_w2=ms_trial(b);
trial_w=unique([trial_w1;trial_w2]);
trial_nw=1:size(tab,1);
trial_nw=setdiff(trial_nw,trial_w)';

acc=tab(:,14);
acc_w=acc(trial_w);
acc_nw=acc(trial_nw);
% acc_w_er=std(acc(trial_w));
% acc_nw_er=std(acc(trial_nw));

rec=tab(:,13)+0.5;
rec_w=rec(trial_w);
rec_nw=rec(trial_nw);

dir=mean(tab(:,12));
tab_w=tab(trial_w,:);
tab_nw=tab(trial_nw,:);
cl_w=tab_w(:,12)>dir;
cl_nw=tab_nw(:,12)>dir;
tab_w_cl=tab_w(cl_w,:);
tab_nw_cl=tab_nw(cl_nw,:);
p_w_hit=sum(tab_w_cl(:,14))/size(tab_w_cl,1);
p_nw_hit=sum(tab_nw_cl(:,14))/size(tab_nw_cl,1);

ncl_w=1:size(tab_w,1);
ncl_w=setdiff(ncl_w,find(cl_w));
ncl_nw=1:size(tab_nw,1);
ncl_nw=setdiff(ncl_nw,find(cl_nw));

tab_w_ncl=tab_w(ncl_w,:);
tab_nw_ncl=tab_nw(ncl_nw,:);
p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);

% d_w=norminv(p_w_hit)-norminv(p_w_fl);
% d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);

%deg=[-8,-4,-2,-1,-0.5,0.5,1,2,4,8];
deg=[0.5,1,2,4,8];
y_w=zeros(10,1);
for i = 1 : 5
    num1=sum(deg(i)==tab_w_ncl(:,11));
    num2=sum(deg(i)==tab_w_cl(:,11));
    y_w(6-i)=num1;
    y_w(5+i)=num2;    
end

y_nw=zeros(10,1);
for i = 1 : 5
    num1=sum(deg(i)==tab_nw_ncl(:,11));
    num2=sum(deg(i)==tab_nw_cl(:,11));
    y_nw(6-i)=num1;
    y_nw(5+i)=num2;    
end
y=[y_w,y_nw];
x=1:10;

figure
bar([1,2],[acc_w,acc_nw])

figure
bar([1,2],[mean(rec_w),mean(rec_nw)])
hold on 
er=errorbar([1,2],[mean(rec_w),mean(rec_nw)],[std(rec_w),std(rec_nw)]);
er.LineStyle= 'none';  
hold off 

% figure
% bar([1,2],[d_w,d_nw])

figure
bar(x,y,'stacked')





end