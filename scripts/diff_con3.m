%radial vs tangenital (original) 
%% this script is changed by plot_figure_diffcon2.m and this is for plot MS count figure
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'UR','LL','UL','LR','VU','VL','HR','HL'};
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
for kk = 1 : 8 
%     Y_con1 = 0 ; 
%     Y_con2 = 0 ; 
%     Y_con3 = 0 ; 
%     Y_con4 = 0 ; 
%     Orig = zeros(10,2);
%     Swit = zeros(10,2);

    tan = zeros(10,2);
    rad = zeros(10,2);
    
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
    
        
    padding_time=[0,100];
    for i = 2 : 2
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
            
            rad_ME = [];
            tan_ME = []; 
            if j ==1 || j == 2 
                for o = 1 : 800
                    if tab(o,9) == 1 || tab(o,9) == 2
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        rad_ME = [rad_ME;m];
                    else
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        tan_ME = [tan_ME;m];  
                    end
                end
            elseif j ==3 || j == 4
                for o = 1 : 800
                    if tab(o,9) == 1 || tab(o,9) == 2
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        tan_ME = [tan_ME;m];
                    else
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        rad_ME = [rad_ME;m];  
                    end
                end                
            elseif j ==5 || j == 6
                for o = 1 : 800
                    if tab(o,9) == 7 || tab(o,9) == 8
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        tan_ME = [tan_ME;m];
                    else
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        rad_ME = [rad_ME;m];  
                    end
                end 
            elseif j ==7 || j == 8
                for o = 1 : 800
                    if tab(o,9) == 5 || tab(o,9) == 6
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        tan_ME = [tan_ME;m];
                    else
                        m=MS_TEMP(MS_TEMP(:,10) == o ,:);
                        rad_ME = [rad_ME;m];  
                    end
                end 
            end
                
                
            [acc_w,acc_nw,rec_w,rec_nw,y_rad,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,rad_ME,samplingRateData,padding_time); 
            [acc_w,acc_nw,rec_w,rec_nw,y_tan,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,tan_ME,samplingRateData,padding_time);             
            tan = tan + y_tan;
            rad = rad + y_rad;          
            
            
    %         tab(:,1)=tab(:,1)+ ((i-1)*8 + (j-1))*800;  
    %         tab_total=[tab_total;tab]

%             [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,padding_time);
%             if j < 5
%                Card = Card + y; 
%             else
%                Obli = Obli + y;
%                 
%             end
            
          
%             ACC_W=[ACC_W;mean(acc_w)];
%             w_num = w_num+size(acc_w,1);
%             ACC_NW=[ACC_NW;mean(acc_nw)];
%             nw_num = nw_num+size(acc_nw,1);

%             REC_W=[REC_W;mean(rec_w(rec_w<4))];
%             REC_NW=[REC_NW;mean(rec_nw(rec_nw<4))];
%              Y=Y+y;
%             p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
%             p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
%             p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
%             p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
% 
%             d_w=norminv(p_w_hit)-norminv(p_w_fl);
%             d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
% 
%             D_w=[D_w;d_w];
%             D_nw=[D_nw;d_nw];
    %         Tab_w_cl=[Tab_w_cl;tab_w_cl];
    %         Tab_nw_cl=[Tab_nw_cl;tab_nw_cl];
    %         Tab_w_ncl=[Tab_w_ncl;tab_w_ncl];
    %         Tab_nw_ncl=[Tab_nw_ncl;tab_nw_ncl];
        end

    
    end
    Orignial=zeros(5,2);
    Switch=zeros(5,2);
    for k = 1 : 5
        Orignial(k,:) = rad(k,:)+rad(11-k,:);
        Switch (k,:)= tan(k,:)+tan(11-k,:);
    end
%     ra = Orignial(:,1)./(Orignial(:,1)+Orignial(:,2));
%     ta = Switch(:,1)./(Switch(:,1)+Switch(:,2));
    ra = Orignial(:,1);
    ta = Switch(:,1);
    Y_total = [ra,ta]/640;
    x=[1,2];
    figure
    c = bar(x,Y_total);
    set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','Original');
    set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','Switched');
    legend(c(1,1:5),'8°','4°','2°','1°','0.5°','Location','eastoutside')
    set(gca,'XTickLabel',{'radial','tangenial'},'FontSize',12);
    ylabel('% MS','FontSize',10);    
    titlename = [subject{kk},' MS rate'];
    title(titlename,'FontSize',12);
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/',subject{kk},'_MS_rate_RvT.png'];
    ylim([0,0.5])
    Original(:,kk) = ra;
    Switched(:,kk) = ta;
    saveas(gca,name)

    
    
    

end

% ori_aab=cumsum(Original,2);
% swi_aab=cumsum(Switched,2);
% x=[0.7,0.85,1,1.15,1.3,1.7,1.85,2,2.15,2.3];
% xt=[1,2];
% Y_total=[ori_aab;swi_aab];
% Oringinal_con = Y_total(1,:) + Y_total(2,:)+ Y_total(3,:)+ Y_total(4,:)+ Y_total(5,:);
% Switch_con = Y_total(6,:) + Y_total(7,:)+ Y_total(8,:)+ Y_total(9,:)+ Y_total(10,:);
% [h,p] = ttest(Oringinal_con,Switch_con);
% output = [Oringinal_con',Switch_con'];
% writematrix(output,'nu_RvT.csv'); 

