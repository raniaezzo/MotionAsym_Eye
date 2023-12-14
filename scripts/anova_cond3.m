clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};

% Y_con1 = 0 ; 
% Y_con2 = 0 ; 
% Y_con3 = 0 ; 
% Y_con4 = 0 ; 
y_stat=[];
y_output=[];
for kk = 1 : 6
    Y_con1 = 0 ; 
    Y_con2 = 0 ; 
    Y_con3 = 0 ; 
    Y_con4 = 0 ; 
    tab_total = [] ;
    MS_total = [] ; 
    sample=[];
%     ACC_W=[];
%     ACC_NW=[];
%     REC_W=[];
%     REC_NW=[];
    ACC_W=zeros(4,4);
    ACC_NW=zeros(4,4);
    REC_W=zeros(4,4);
    REC_NW=zeros(4,4);
    Y=zeros(10,2);
    % Tab_w_cl=[];
    % Tab_nw_cl=[];
    % Tab_w_ncl=[];
    % Tab_nw_ncl=[];
    D_w=zeros(4,4);
    D_nw=zeros(4,4);
   %calculate numbers of different condition
    w_num_con1 = 0;
    nw_num_con1= 0;
    w_num_con2 = 0;
    nw_num_con2= 0;
    w_num_con3 = 0;
    nw_num_con3= 0;
    w_num_con4 = 0;
    nw_num_con4= 0;

%     Y_con1 = 0 ; 
%     Y_con2 = 0 ; 
%     Y_con3 = 0 ; 
%     Y_con4 = 0 ; 
        
    padding_time=[0,100];
    for i = 1 : 2
        for j = 1 : 8

            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{kk}, ...
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
    %         tab_total=[tab_total;tab]

            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,padding_time);
            y1 = (y([1,2,3,4,5],1)+y([10,9,8,7,6],1))/160;
            y_stat=[y_stat;y1];
          

        end
    end

end

%% Sat
sub_stat = [ones(80,1);ones(80,1)*2;ones(80,1)*3;ones(80,1)*4;ones(80,1)*5;ones(80,1)*6];
fa1 = [8;4;2;1;0.5];
fa1_stat = [];
for aa = 1 :96
   fa1_stat = [fa1_stat;fa1];
end
polar = [ones(40,1)*2;ones(40,1)];% 1 for original, 2 for switched
polar_stat=[];
for bb = 1 :6
   polar_stat = [polar_stat;polar];
end


Car = [ones(20,1);ones(20,1)*2]; % 1 for cardinal, 2 for oblique
car_stat=[];
for cc = 1 :12
   car_stat = [car_stat;Car];
end

y_stat=(y_stat - min(y_stat))./(max(y_stat)-min(y_stat)); %normalize to [0,1]

mat = [y_stat,sub_stat,fa1_stat,polar_stat,car_stat];

%% Overall reult
reg = fitlm(mat(:,3),mat(:,1));
figure
plot(reg)
xlabel('Difficult')
ylabel('MS rate')
title('All Subject')
figpath='C:\Users\86186\Desktop\fig';
name = [figpath,'/','_All_reg.png'];
saveas(gca,name)

names = {'_sub1_reg.png','_sub2_reg.png','_sub3_reg.png','_sub4_reg.png','_sub5_reg.png','_sub6_reg.png'};

for as = 1 : 6
    m=mat(mat(:,2)==as,:);
    re = fitlm(m(:,3),m(:,1));
    figure
    plot(re)
    xlabel('Difficult')
    ylabel('MS rate')
    title(subject{as})
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/',names{as}];
    saveas(gca,name)
end

%% cardinal vs oblique (Polar)
rad = mat(mat(:,4)==1,:);
tan = mat(mat(:,4)==2,:);
rad_all = fitlm(rad(:,3),rad(:,1));
tan_all = fitlm(tan(:,3),tan(:,1));

% figure
% %subplot(2,1,1)
% plot(rad_all)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Cardinal(Polar)')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','cardinal_polar_reg.png'];
% saveas(gca,name)
% 
% figure
% %subplot(2,1,2)
% plot(tan_all)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Oblique(Polar)')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','oblique_polar_reg.png'];
% saveas(gca,name)
rad_val_x1P = zeros(6,1);
tan_val_x1P = zeros(6,1);
rad_val_intP = zeros(6,1);
tan_val_intP = zeros(6,1);
for i = 1 : 6 
    r= rad(rad(:,2)==i,:);
    t= tan(tan(:,2)==i,:);
    rad_all = fitlm(r(:,3),r(:,1));
    tan_all = fitlm(t(:,3),t(:,1));
    cdc_rad = rad_all.Coefficients;
    cdc_tan = tan_all.Coefficients;
    rad_val_x1P(i) = cdc_rad.Estimate(2);
    tan_val_x1P(i) = cdc_tan.Estimate(2);
    rad_val_intP(i) = cdc_rad.Estimate(1);
    tan_val_intP(i) = cdc_tan.Estimate(1);
%     figure
%     %subplot(2,1,1)
%     plot(rad_all)
%     xlabel('Difficult')
%     ylabel('MS rate')
%     titlename = ['Cardinal Sub',num2str(i)];
%     title(titlename)
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/','Cardinal_Polar',num2str(i),'_reg.png'];
%     saveas(gca,name)
% 
%     figure
%     %subplot(2,1,2)
%     plot(tan_all)
%     xlabel('Difficult')
%     ylabel('MS rate')
%     titlename = ['Oblique Sub',num2str(i)];
%     title(titlename)
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/','Oblique_Polar',num2str(i),'_reg.png'];
%     saveas(gca,name)
end
[h_radv0P,p_radv0P] = ttest(rad_val_x1P);
[h_tanv0P,p_tanv0P] = ttest(tan_val_x1P);
[h_tanvradP,p_tanvradP] = ttest(tan_val_x1P,rad_val_x1P);
output_PCVCO = [rad_val_x1P,tan_val_x1P];
% writematrix(output_PCVCO,'OvS.csv'); 
%% cardinal vs oblique (Cartesian)
rad = mat(mat(:,5)==1,:);
tan = mat(mat(:,5)==2,:);
rad_all = fitlm(rad(:,3),rad(:,1));
tan_all = fitlm(tan(:,3),tan(:,1));

% figure
% %subplot(2,1,1)
% plot(rad_all)
% xlabel('Difficult')
% ylabel('MS rate')
% title('cardinal')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','cardinal(Cartesian).png'];
% saveas(gca,name)
% 
% figure
% %subplot(2,1,2)
% plot(tan_all)
% xlabel('Difficult')
% ylabel('MS rate')
% title('oblique')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','oblique(Cartesian).png'];
% saveas(gca,name)
rad_val_x1 = zeros(6,1);
tan_val_x1 = zeros(6,1);
rad_val_int = zeros(6,1);
tan_val_int = zeros(6,1);

for i = 1 : 6 
    r= rad(rad(:,2)==i,:);
    t= tan(tan(:,2)==i,:);
    rad_all = fitlm(r(:,3),r(:,1));
    tan_all = fitlm(t(:,3),t(:,1));
    cdc_rad = rad_all.Coefficients;
    cdc_tan = tan_all.Coefficients;
    rad_val_x1(i) = cdc_rad.Estimate(2);
    tan_val_x1(i) = cdc_tan.Estimate(2);
    rad_val_int(i) = cdc_rad.Estimate(1);
    tan_val_int(i) = cdc_tan.Estimate(1);
    figure
    %subplot(2,1,1)
    plot(rad_all)
    xlabel('Difficult')
    ylabel('MS rate')
    titlename = ['Cardinal Sub',num2str(i)];
    title(titlename)
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/','Cardinal_Cartesian',num2str(i),'_reg.png'];
    saveas(gca,name)

    figure
    %subplot(2,1,2)
    plot(tan_all)
    xlabel('Difficult')
    ylabel('MS rate')
    titlename = ['Oblique Sub',num2str(i)];
    title(titlename)
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/','Oblique_Cartesian',num2str(i),'_reg.png'];
    saveas(gca,name)
end

[h_radv0,p_radv0] = ttest(rad_val_x1);
[h_tanv0,p_tanv0] = ttest(tan_val_x1);
[h_tanvrad,p_tanvrad] = ttest(tan_val_x1,rad_val_x1);

output_CCVCO = [rad_val_x1,tan_val_x1];
% writematrix(output_CCVCO,'CvT.csv');
