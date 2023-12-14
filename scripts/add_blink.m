%% this function add blink into tab file
function path=add_blink(tab_path,blink_path)
[path,filename, ~] = fileparts(tab_path);
load(tab_path)
load(blink_path)
tab_len=size(tab,1);
c1=zeros(tab_len,1);
for i = 1 : tab_len
    order=tab(i,2)+1300 <=blink(:,1) & tab(i,2)+1800 >=blink(:,1);
    c1(i)=sum(order);
end
tab=[tab,c1];
save(sprintf('%s/%s_blink.mat',path,filename),'tab');
path=sprintf('%s/%s_blink.mat',path,filename);
end