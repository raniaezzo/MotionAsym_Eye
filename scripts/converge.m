%% create a new tab file, which contain behavior data
%filepath1= tab filepath way/ filepath2= csv filepath way/
%filepath3 = Dat_stim file path
function path=converge(filepath1,filepath2)
tab_stim=load(filepath1);
M_csv=csvread(filepath2);
%Dat_stim=load(filepath3);
tab=tab_stim.tab;
a=tab(:,7)==0; 
tab=tab(a,:);
M_csv=M_csv(a,:);
c=M_csv(:,3);
d=M_csv(:,4);
f=M_csv(:,6);
j=M_csv(:,10);
m=M_csv(:,13);
n=M_csv(:,14);
tab_new=[tab c d f j m n];
[path, filename, ~] = fileparts(filepath1);
save(sprintf('%s/%s_new.mat',path,filename),'tab_new');
path=sprintf('%s/%s_new.mat',path,filename);
end