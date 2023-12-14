% Script for first sperate with and not with, and seperate tilt angle and
% calculate accuracy 

clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
y_cond1=zeros(5,6);
y_cond2=zeros(5,6);
y_cond3=zeros(5,6);
y_cond4=zeros(5,6);
Y_con1 = 0 ; 
Y_con2 = 0 ; 
Y_con3 = 0 ; 
Y_con4 = 0 ; 
Original = zeros(5,6);
Switched = zeros(5,6);
tilt = [8,4,2,1,0.5];
acc_with_sub =zeros(5,6);
acc_notwith_sub =zeros(5,6);
ACC = zeros(2,6);
for kk = 1 : 6
    tab_nw_total = [];
    tab_w_total = [];
%     Y_con1 = 0 ; 
%     Y_con2 = 0 ; 
%     Y_con3 = 0 ; 
%     Y_con4 = 0 ; 
    Orig = zeros(10,2);
    Swit = zeros(10,2);
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
    acc_with =zeros(5,1);
    acc_notwith = zeros(5,1);
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
            
            tab_nw = [tab_nw_cl;tab_nw_ncl];
            tab_w = [tab_w_cl;tab_w_ncl];
            tab_nw_total = [tab_nw_total;tab_nw];
            tab_w_total = [tab_w_total;tab_w];
            
            
                
          

        end

    
    end
    acc = zeros(1,2);
    acc(1) = mean(tab_w_total(:,15));
    acc(2) = mean(tab_nw_total(:,15));
    ACC(:,kk) = acc';
%     x = 1:2;
%     figure
%     bar(x,acc)
%     ylim([0.5,1.1])
% %    legend('With MS','Not with MS')
%     xticklabels({'With MS','Not with MS'})
%     ylabel('Accuracy')
%     titlename = [subject{kk},' Accuracy'];
%     title(titlename,'FontSize',12);
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/',subject{kk},'_sumacc.png'];
%     saveas(gca,name)


%     for ti = 1:5
%         a_nw = mean(tab_nw_total(tab_nw_total(:,11)==tilt(ti),15));
%         a_w = mean(tab_w_total(tab_w_total(:,11)==tilt(ti),15));
%         acc_with(ti)=a_w;
%         acc_notwith(ti) = a_nw;
%     end
%     a_wnw = [acc_with,acc_notwith]';
%     x = 1:5;
%     figure
%     bar(x,a_wnw)
%     ylim([0.5,1.1])
%     legend('With MS','Not with MS')
%     xticklabels({'8°','4°','2°','1°','0.5°'})
%     ylabel('Accuracy')
%     titlename = [subject{kk},' Accuracy'];
%     title(titlename,'FontSize',12);
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/',subject{kk},'_totalacc.png'];
%     saveas(gca,name)
%     
%     acc_with_sub(:,kk)=acc_with;
%     acc_notwith_sub(:,kk)=acc_notwith;
    

%     Orignial=zeros(5,2);
%     Switch=zeros(5,2);
%     for k = 1 : 5
%         Orignial(k,:) = Orig(k,:)+Orig(11-k,:);
%         Switch (k,:)= Swit(k,:)+Swit(11-k,:);
%     end
%     or = Orignial(:,1)./(Orignial(:,1)+Orignial(:,2));
%     sw = Switch(:,1)./(Switch(:,1)+Switch(:,2));
%     or = Orignial(:,1);
%     sw = Switch(:,1);
%     Y_total = [or,sw]/1280;
%     x=[1,2];
%     figure
%     c = bar(x,Y_total);
%     set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','Original');
%     set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','Switched');
%     lgd=legend(c(1,1:5),'8°','4°','2°','1°','0.5°','Location','eastoutside');
% 
%     set(gca,'XTickLabel',{'Polar Cardinal','Polar Oblique'},'FontSize',12);
%     ylabel('% MS','FontSize',10);    
%     titlename = [subject{kk},' MS rate'];
%     title(titlename,'FontSize',12);
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/',subject{kk},'_MS_rate_OvS.png'];
%     ylim([0,0.5])
%     Original(:,kk) = or;
%     Switched(:,kk) = sw;
%     saveas(gca,name)

    
     
    
    

end
x = 1 :6;
% A = zeros(2,5);
% A(1,:)=mean(acc_with_sub,2)';
% A(2,:)=mean(acc_notwith_sub,2)';
figure
bar(x,ACC)
ylim([0.5,1.0])
legend('With MS','Not with MS')
xticklabels({'sub1','sub2','sub3','sub4','sub5','sub6'})
ylabel('Accuracy')
titlename = ['Accuracy'];
title(titlename,'FontSize',12);
figpath='C:\Users\86186\Desktop\fig';
name = [figpath,'/','subtilt.png'];
% saveas(gca,name)