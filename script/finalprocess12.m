clc; clear all;

tab_total = [] ;
MS_total = [] ; 
sample=[];
sub={'S01','S02','S03','S04','S05','S06','S07','S08'};
subject = 'S01'; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};

Y=zeros(10,2);


acc_total_w=[];
acc_total_nw=[];
rec_total_w=[];
rec_total_nw=[];
d_total_w=[];
d_total_nw=[];
figpath='C:\Users\86186\Desktop\fig\new';

color={[215,25,28],[253,174,97],[255,255,191],[171,217,233],[44,123,182]};	
y_total = zeros(8,10);
for k = 1 : 8
    Y=zeros(10,2);
%    y_total=0;
    for i = 1 : 2
        if i ==1 & k > 6
            continue
        end
        for j = 1 : 8

            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', sub{k}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
            sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
    %         cd(fullfile(main_folder, 'eyedata','MATs'));
    %         ms_name = dir(sprintf('%s.mat', direction{j})).name;
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction{j}) );
            load(ms_path);
    %         MS_TEMP(:,10)=MS_TEMP(:,10)+ ((i-1)*8 + (j-1))*800;   
    %         MS_total=[MS_total;MS_TEMP];
    %         MATpath = fullfile(main_folder, 'eyedata','MATs');
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
    %         tab(:,1)=tab(:,1)+ ((i-1)*8 + (j-1))*800;  
    %         tab_total=[tab_total;tab];
            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,[0,100]);

            Y=Y+y;

        end


    end

    y_total(k,:) = Y(:,1)';

    
end


y_plot = zeros(8,5);
for sd = 1 :5 
    y_plot(:,sd) = y_total(:,sd) + y_total(:,(11-sd));
end

y_plot_nor = zeros(8,5);
y_plot_nor(1:6,:) = y_plot(1:6,:)/2560;
y_plot_nor([7,8],:) = y_plot([7,8],:)/1280;

y_nor = mean(y_plot_nor,1);


X = categorical({'8','4','2','1','0.5'});
X = reordercats(X,{'8','4','2','1','0.5'});
figure
b = bar(X,y_nor);
b.FaceColor = 'flat';
for sdf = 1 : 5 
    b.CData(sdf,:) = color{6-sdf}/255;
end
xlabel('tilt angle')
ylabel('% of trial')
title('% of trial with microsaccade during stimulus')
name = 'C:\Users\86186\Desktop\fig\new\tilt_nor.png';
saveas(gca,name)