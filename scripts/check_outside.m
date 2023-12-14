%% check the outside points 
%this function can save the clear data and output the number of points for
%each trial, and the input is the filepath for new Dat_stim data
function path=check_outside(filepath)
[path, filename,~] = fileparts(filepath);
load(filepath);
aa=size(Dat_stim,1);
outside=zeros(aa,2);
for jj = 1 : aa 
    outside(jj,1)=(Dat_stim(jj,2)-576)^2+(Dat_stim(jj,3)-435)^2 > 59.33^2;
end
outside(:,2)=Dat_stim(:,1);
save(sprintf('%s/%s_outside.mat',path,filename),'outside')
path=sprintf('%s/%s_outside.mat',path,filename);
end