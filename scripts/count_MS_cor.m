%% this function is used to y
function y= count_MS_cor(tab,MS_TEMP,samplerate,time)
lower=(1300-time(1))*samplerate/1000;
higher=(1800+time(2))*samplerate/1000;
ms_start=MS_TEMP(:,8);
ms_end=MS_TEMP(:,9);
ms_trial=MS_TEMP(:,10);
a= lower<ms_start & higher > ms_start;
trial_w1=ms_trial(a);
b= lower<ms_end & higher > ms_end;
trial_w2=ms_trial(b);
trial_w=unique([trial_w1;trial_w2]); % w MS
trial_nw=1:size(tab,1);
trial_nw=setdiff(trial_nw,trial_w)'; % nw MS
tab_w=[];
tab_nw=[];
for aa = trial_w'
    tab_w=[tab_w;tab(tab(:,1)==aa,:)];
end
for bb = trial_nw'
    tab_nw=[tab_nw;tab(tab(:,1)==bb,:)];
end


cl_w=tab_w(:,12)==1; % clockwise w MS
cl_nw=tab_nw(:,12)==1; % clockwise w/o MS (corrected BUG from 0 to 1)
tab_w_cl=tab_w(cl_w,:); % clockwise w MS
tab_nw_cl=tab_nw(cl_nw,:); % counterclockwise w/o MS

p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1); % hits with MS
p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1); % hits w/o MS

ncl_w=1:size(tab_w,1); 
ncl_w=setdiff(ncl_w,find(cl_w)); % counterclock w MS
ncl_nw=1:size(tab_nw,1); 
ncl_nw=setdiff(ncl_nw,find(cl_nw)); % counterclock w/o MS

tab_w_ncl=tab_w(ncl_w,:);
tab_nw_ncl=tab_nw(ncl_nw,:);
p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);

d_w=norminv(p_w_hit)-norminv(p_w_fl);
d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);

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
end
