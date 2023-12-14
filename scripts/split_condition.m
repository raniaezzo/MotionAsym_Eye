%% split each condition
%we input the new tab file and Data_stim filepath and we save each
%condition into new data file. filepath1---tab file, filepath2---Data_stim
%file
function split_condition(filepath1,filepath2)
loc={'loc_315','loc_135','loc_225','loc_45','loc_270','loc_90','loc_180','loc_0'};
% [path1, filename1,~] = fileparts(filepath1);
[path2, filename2,~] = fileparts(filepath2);
file='seperate condition';
Dat_stim=load(filepath2);
tab_stim=load(filepath1);
Dat_stim=Dat_stim.Dat_stim;
tab=tab_stim.tab_new;
locat_num=unique(tab(:,9));
tab_con=zeros(200,14,4);
for f = 1:4 %4 locations
    b = tab(:,9) == locat_num(f);
    tab_con(:,:,f)=tab(b,:);
end
num_timepoints=sum(Dat_stim(:,1)>=tab(1,2,1)+1300 & Dat_stim(:,1)<=tab(1,2,1)+1800);
final_x = nan(4,200,num_timepoints+1);
final_y = nan(4,200,num_timepoints+1);

for j = 1 : 4 
    for k = 1 : 200
        num=tab_con(k,5,j);
        order=num<=Dat_stim(:,1) & Dat_stim(:,1)<=num+500;        
        x=Dat_stim(order,2);
        y=Dat_stim(order,3);
        final_x(j,k,1:sum(order))=squeeze(squeeze(x));
        final_y(j,k,1:sum(order))=squeeze(squeeze(y)); 
     end
end

for ii = 1 : 4 
    if locat_num(1) == 1
        con=squeeze(final_x(ii,:,:));
        name=[loc{ii},filename2];
        save(sprintf('%s/%s/%s_new.mat',path2,file,name),'con')
    elseif locat_num(1) == 5
        con=squeeze(final_x(ii,:,:));
        name=[loc{ii+4},'_',filename2];
        save(sprintf('%s/%s/%s_new.mat',path2,file,name),'con')
    end
end
end