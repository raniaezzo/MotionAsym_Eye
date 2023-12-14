%histogram figure
function draw_his(ms_new,filepath)
%bound=-157.5:45:157.5;
sum_num=length(ms_new(:,13));
[path,filename, ~] = fileparts(filepath);
edge=0.125:0.25:2.125;
edge=edge-1;
edge1=edge*pi;

num_com=zeros(1,8);
num_am=zeros(1,8);
for i = 1 : 7 
    num_com(i)=sum(bound(i)<ms_new(:,13) & bound(i+1)>ms_new(:,13));
    num_am(i)=sum(bound(i)<ms_new(:,12) & bound(i+1)>ms_new(:,12));
end
num_com(8)=sum_num-sum(num_com);
num_am(8)=sum_num-sum(num_am);
a=figure
subplot(1,2,1)
polarhistogram('BinEdges',edge1,'BinCounts',num_com)
title('component')
subplot(1,2,2)
%polarhistogram('BinEdges',edge1,'BinCounts',num_am)
polarhistogram(deg2rad(ms_new(:,13)),deg2rad([-22.5 22.5 67.5 112.5 157.5 202.5 247.5 292.5 337.49999])),
title('amplitude')

savefig(a,sprintf('%s/%s_figure3.fig',path,filename))
end