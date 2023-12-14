%% this function to draw a different figure
function countfigure(MS,tab,samplerate,filepath)
locationids = 1:8; locationdegrees = {315,135,225,45,270,90,180,0};
locationlabels = strcat('loc',cellfun(@num2str,locationdegrees,'un',0));
[path,filename, ~] = fileparts(filepath);
%%should add samplerate
x=(1000/samplerate):(1000/samplerate):2500;
count=size(MS,1);
timepoint=zeros(count,2500*samplerate/1000);
for i = 1 : count
    a=MS(i,8);
    b=MS(i,9);
    timepoint(i,a:b)=ones(1,b-a+1);
end
countnum=sum(timepoint,1);

c=zeros(count,1);
for j = 1 : 800
    con=tab(j,9);
    if con == 1 ||  con == 5
        aa=MS(:,10)==j;
        c(aa)=1;
    elseif con == 2 ||  con == 6
        aa=MS(:,10)==j;
        c(aa)=2;
    elseif con == 3  ||  con == 7 
        aa=MS(:,10)==j;
        c(aa)=3;
    else
        aa=MS(:,10)==j;
        c(aa)=4;
    end
end

c1=sum(timepoint(c==1,:),1);
c2=sum(timepoint(c==2,:),1);
c3=sum(timepoint(c==3,:),1);
c4=sum(timepoint(c==4,:),1);

a=figure
subplot(2,1,1)
plot(x,countnum)
subplot(2,1,2)
hold on 
plot(x,c1)
hold on 
plot(x,c2)
hold on 
plot(x,c3)
hold on 
plot(x,c4)
hold on 
loca=tab(1,9);
if loca == 1 || loca == 2 || loca == 3 || loca == 4
    legend(locationlabels{1},locationlabels{2},locationlabels{3},locationlabels{4})
else
    legend(locationlabels{5},locationlabels{6},locationlabels{7},locationlabels{8})
end
savefig(a,sprintf('%s/%s_figure2.fig',path,filename))   
end