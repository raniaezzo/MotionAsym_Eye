%% copy from fig_new, and it is used to plot the correct and incorrect figure
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
y_cond1=zeros(5,6);
y_cond2=zeros(5,6);
y_cond3=zeros(5,6);
y_cond4=zeros(5,6);
% Y_con1_co = 0 ; 
% Y_con2_co = 0 ; 
% Y_con3_co = 0 ; 
% Y_con4_co = 0 ; 
% Y_con1_unco= 0 ; 
% Y_con2_unco = 0 ; 
% Y_con3_unco = 0 ; 
% Y_con4_unco = 0 ; 
for kk = 1 : 6 
    Y_con1_co = 0 ; 
    Y_con2_co = 0 ; 
    Y_con3_co = 0 ; 
    Y_con4_co = 0 ; 
    Y_con1_unco= 0 ; 
    Y_con2_unco = 0 ; 
    Y_con3_unco = 0 ; 
    Y_con4_unco = 0 ; 
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

%            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,padding_time);
            %% split correct and incorrect
            tab_cor = tab(tab(:,15)==1,:);
            tab_incor = tab(tab(:,15)==0,:);
            MS_cor = [];
            MS_incor = [];
            ord_cor = find(tab(:,15)==1);
            ord_incor=find(tab(:,15)==0);
            for cor = ord_cor'
                MS_cor = [MS_cor;MS_TEMP(MS_TEMP(:,10)==cor,:)];
            end
            for incor=ord_incor'
                MS_incor = [MS_incor;MS_TEMP(MS_TEMP(:,10)==incor,:)];
            end
            y_cor= count_MS_cor(tab_cor,MS_cor,samplingRateData,padding_time);
            y_incor= count_MS_cor(tab_incor,MS_incor,samplingRateData,padding_time);
            %%
            if i == 2 
                if j < 5
                    Y_con1_co = Y_con1_co+y_cor;
                    Y_con1_unco = Y_con1_unco+y_incor;
                else
                 
                    Y_con2_co = Y_con2_co+y_cor;
                    Y_con2_unco = Y_con2_unco+y_incor;
                end
            else
                if j < 5
                   
                    Y_con3_co = Y_con3_co+y_cor;
                    Y_con3_unco = Y_con3_unco+y_incor;
                else
                    
                    Y_con4_co = Y_con4_co+y_cor;
                    Y_con4_unco = Y_con4_unco+y_incor;
                end
            end
          

        end

    
    end
%% for each subject's figure
    y_con1_cor = Y_con1_co([1,2,3,4,5],:)+Y_con1_co([10,9,8,7,6],:);
    y_con2_cor = Y_con2_co([1,2,3,4,5],:)+Y_con2_co([10,9,8,7,6],:);
    y_con3_cor = Y_con3_co([1,2,3,4,5],:)+Y_con3_co([10,9,8,7,6],:);
    y_con4_cor = Y_con4_co([1,2,3,4,5],:)+Y_con4_co([10,9,8,7,6],:);
    y_con1_uncor = Y_con1_unco([1,2,3,4,5],:)+Y_con1_unco([10,9,8,7,6],:);
    y_con2_uncor = Y_con2_unco([1,2,3,4,5],:)+Y_con2_unco([10,9,8,7,6],:);
    y_con3_uncor = Y_con3_unco([1,2,3,4,5],:)+Y_con3_unco([10,9,8,7,6],:);
    y_con4_uncor = Y_con4_unco([1,2,3,4,5],:)+Y_con4_unco([10,9,8,7,6],:);
    %calculate rate
    r_con1_cor=y_con1_cor(:,1)./(y_con1_cor(:,1)+y_con1_cor(:,2));
    r_con2_cor=y_con2_cor(:,1)./(y_con2_cor(:,1)+y_con2_cor(:,2));
    r_con3_cor=y_con3_cor(:,1)./(y_con3_cor(:,1)+y_con3_cor(:,2));
    r_con4_cor=y_con4_cor(:,1)./(y_con4_cor(:,1)+y_con4_cor(:,2));

    r_con1_uncor=y_con1_uncor(:,1)./(y_con1_uncor(:,1)+y_con1_uncor(:,2));
    r_con2_uncor=y_con2_uncor(:,1)./(y_con2_uncor(:,1)+y_con2_uncor(:,2));
    r_con3_uncor=y_con3_uncor(:,1)./(y_con3_uncor(:,1)+y_con3_uncor(:,2));
    r_con4_uncor=y_con4_uncor(:,1)./(y_con4_uncor(:,1)+y_con4_uncor(:,2));

    figure
    x=[1 2 3 4];
    r_con_cor=[r_con1_cor,r_con2_cor,r_con3_cor,r_con4_cor];
    c=bar(x,r_con_cor');
    set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','8');
    set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','4');
    set(c(1,3),'FaceColor','b','BarWidth',0.9,'DisplayName','2');
    set(c(1,4),'FaceColor','r','BarWidth',0.9,'DisplayName','1');
    set(c(1,5),'FaceColor','g','BarWidth',0.9,'DisplayName','0.5');
    set(gca,'XTickLabel',{'Cardinal directions ON meridians','Oblique directions OFF meridians','Cardinal directions OFF meridians','Oblique directions ON meridians'},'FontSize',7);
    legend(c(1,1:5),'8°','4°','2°','1°','0.5°');
    ylabel('% trial','FontSize',10);
    ylim([0,0.5]);
    titlename = [subject{kk},' MS rate (correct trial)'];
    title(titlename,'FontSize',12);
    figpath='C:\Users\86186\Desktop\fig code';
    name = [figpath,'/',subject{kk},'_MS_rate_diffcon_corr.png'];
    saveas(gca,name)

    figure
    x=[1 2 3 4];
    r_con_uncor=[r_con1_uncor,r_con2_uncor,r_con3_uncor,r_con4_uncor];
    c=bar(x,r_con_uncor');
    set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','8');
    set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','4');
    set(c(1,3),'FaceColor','b','BarWidth',0.9,'DisplayName','2');
    set(c(1,4),'FaceColor','r','BarWidth',0.9,'DisplayName','1');
    set(c(1,5),'FaceColor','g','BarWidth',0.9,'DisplayName','0.5');
    set(gca,'XTickLabel',{'Cardinal directions ON meridians','Oblique directions OFF meridians','Cardinal directions OFF meridians','Oblique directions ON meridians'},'FontSize',7);
    legend(c(1,1:5),'8°','4°','2°','1°','0.5°');
    ylabel('% trial','FontSize',10);
    ylim([0,0.5]);
    titlename = [subject{kk},' MS rate (incorrect trial)'];
    title(titlename,'FontSize',12);
    figpath='C:\Users\86186\Desktop\fig code';
    name = [figpath,'/',subject{kk},'_MS_rate_diffcon_incorr.png'];
    saveas(gca,name)
end
%% for all subjects figure
%combine + and - 
% y_con1_cor = Y_con1_co([1,2,3,4,5],:)+Y_con1_co([10,9,8,7,6],:);
% y_con2_cor = Y_con2_co([1,2,3,4,5],:)+Y_con2_co([10,9,8,7,6],:);
% y_con3_cor = Y_con3_co([1,2,3,4,5],:)+Y_con3_co([10,9,8,7,6],:);
% y_con4_cor = Y_con4_co([1,2,3,4,5],:)+Y_con4_co([10,9,8,7,6],:);
% y_con1_uncor = Y_con1_unco([1,2,3,4,5],:)+Y_con1_unco([10,9,8,7,6],:);
% y_con2_uncor = Y_con2_unco([1,2,3,4,5],:)+Y_con2_unco([10,9,8,7,6],:);
% y_con3_uncor = Y_con3_unco([1,2,3,4,5],:)+Y_con3_unco([10,9,8,7,6],:);
% y_con4_uncor = Y_con4_unco([1,2,3,4,5],:)+Y_con4_unco([10,9,8,7,6],:);
% %calculate rate
% r_con1_cor=y_con1_cor(:,1)./(y_con1_cor(:,1)+y_con1_cor(:,2));
% r_con2_cor=y_con2_cor(:,1)./(y_con2_cor(:,1)+y_con2_cor(:,2));
% r_con3_cor=y_con3_cor(:,1)./(y_con3_cor(:,1)+y_con3_cor(:,2));
% r_con4_cor=y_con4_cor(:,1)./(y_con4_cor(:,1)+y_con4_cor(:,2));
% 
% r_con1_uncor=y_con1_uncor(:,1)./(y_con1_uncor(:,1)+y_con1_uncor(:,2));
% r_con2_uncor=y_con2_uncor(:,1)./(y_con2_uncor(:,1)+y_con2_uncor(:,2));
% r_con3_uncor=y_con3_uncor(:,1)./(y_con3_uncor(:,1)+y_con3_uncor(:,2));
% r_con4_uncor=y_con4_uncor(:,1)./(y_con4_uncor(:,1)+y_con4_uncor(:,2));
% 
% figure
% x=[1 2 3 4];
% r_con_cor=[r_con1_cor,r_con2_cor,r_con3_cor,r_con4_cor];
% c=bar(x,r_con_cor');
% set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','8');
% set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','4');
% set(c(1,3),'FaceColor','b','BarWidth',0.9,'DisplayName','2');
% set(c(1,4),'FaceColor','r','BarWidth',0.9,'DisplayName','1');
% set(c(1,5),'FaceColor','g','BarWidth',0.9,'DisplayName','0.5');
% set(gca,'XTickLabel',{'Cardinal directions ON meridians','Oblique directions OFF meridians','Cardinal directions OFF meridians','Oblique directions ON meridians'},'FontSize',7);
% legend(c(1,1:5),'8°','4°','2°','1°','0.5°');
% ylabel('% trial','FontSize',10);
% ylim([0,0.5])
% title('MS rate (correct trial)','FontSize',12);
% figpath='C:\Users\86186\Desktop\fig code';
% name = [figpath,'/','MS_rate_diffcon_corr.png'];
% saveas(gca,name)
% 
% figure
% x=[1 2 3 4];
% r_con_uncor=[r_con1_uncor,r_con2_uncor,r_con3_uncor,r_con4_uncor];
% c=bar(x,r_con_uncor');
% set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','8');
% set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','4');
% set(c(1,3),'FaceColor','b','BarWidth',0.9,'DisplayName','2');
% set(c(1,4),'FaceColor','r','BarWidth',0.9,'DisplayName','1');
% set(c(1,5),'FaceColor','g','BarWidth',0.9,'DisplayName','0.5');
% set(gca,'XTickLabel',{'Cardinal directions ON meridians','Oblique directions OFF meridians','Cardinal directions OFF meridians','Oblique directions ON meridians'},'FontSize',7);
% legend(c(1,1:5),'8°','4°','2°','1°','0.5°');
% ylabel('% trial','FontSize',10);
% ylim([0,0.5]);
% title('MS rate (incorrect trial)','FontSize',12);
% figpath='C:\Users\86186\Desktop\fig code';
% name = [figpath,'/','MS_rate_diffcon_incorr.png'];
% saveas(gca,name)



