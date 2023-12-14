clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'UR','LL','UL','LR','VU','VL','HR','HL'};

% Y_con1 = 0 ; 
% Y_con2 = 0 ; 
% Y_con3 = 0 ; 
% Y_con4 = 0 ; 
y_stat=[];
y_output=[];
T_y = [];
for kk = 1 : 8
    tan = zeros(10,2);
    rad = zeros(10,2);
    Y_rad_con1 = 0 ; 
    Y_rad_con2 = 0 ; 
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
    %         tab(:,1)=tab(:,1)+ ((i-1)*8 + (j-1))*800;  
    %         tab_total=[tab_total;tab]
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
%[acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,padding_time);
%             if i == 2 
%                 if j < 5
%                     %acc
% %                     ACC_W(1,j)=mean(acc_w);
% %                     ACC_NW(1,j)=mean(acc_nw);
% %                     w_num_con1 = w_num_con1+size(acc_w,1);
% %                     nw_num_con1= nw_num_con1+size(acc_nw,1);
% %                     %rec
% %                     REC_W(1,j)=mean(rec_w(rec_w<4));
% %                     REC_NW(1,j)=mean(rec_nw(rec_nw<4));
% %                     %d_prime
% %                     p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
% %                     p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
% %                     p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
% %                     p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
% %                     d_w=norminv(p_w_hit)-norminv(p_w_fl);
% %                     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
% %                     D_w(1,j)=d_w;
% %                     D_nw(1,j)=d_nw;
%                     Y_rad_con1 = Y_rad_con1+y_rad;
%                     Y_tan_con1 = Y_tan_con1+y_tan;
%                 else
% %                     ACC_W(2,j-4)=mean(acc_w);
% %                     ACC_NW(2,j-4)=mean(acc_nw);
% %                     w_num_con2 = w_num_con2+size(acc_w,1);
% %                     nw_num_con2= nw_num_con2+size(acc_nw,1);
% %                     %rec
% %                     REC_W(2,j-4)=mean(rec_w(rec_w<4));
% %                     REC_NW(2,j-4)=mean(rec_nw(rec_nw<4));
% %                     %d_prime
% %                     p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
% %                     p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
% %                     p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
% %                     p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
% %                     d_w=norminv(p_w_hit)-norminv(p_w_fl);
% %                     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
% %                     D_w(2,j-4)=d_w;
% %                     D_nw(2,j-4)=d_nw;
%                     Y_rad_con2 = Y_rad_con2+y_rad;
%                     Y_tan_con2 = Y_tan_con1+y_tan;
%                 end
%             else
%                 if j < 5
%                     ACC_W(3,j)=mean(acc_w);
%                     ACC_NW(3,j)=mean(acc_nw);
%                     w_num_con3 = w_num_con3+size(acc_w,1);
%                     nw_num_con3= nw_num_con3+size(acc_nw,1);
%                     %rec
%                     REC_W(3,j)=mean(rec_w(rec_w<4));
%                     REC_NW(3,j)=mean(rec_nw(rec_nw<4));
%                     %d_prime
%                     p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
%                     p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
%                     p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
%                     p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
%                     d_w=norminv(p_w_hit)-norminv(p_w_fl);
%                     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
%                     D_w(3,j)=d_w;
%                     D_nw(3,j)=d_nw;
%                     Y_con3 = Y_con3+y;
%                 else
%                     ACC_W(4,j-4)=mean(acc_w);
%                     ACC_NW(4,j-4)=mean(acc_nw);
%                     w_num_con4 = w_num_con4+size(acc_w,1);
%                     nw_num_con4= nw_num_con4+size(acc_nw,1);
%                     %rec
%                     REC_W(4,j-4)=mean(rec_w(rec_w<4));
%                     REC_NW(4,j-4)=mean(rec_nw(rec_nw<4));
%                     %d_prime
%                     p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
%                     p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
%                     p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
%                     p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
%                     d_w=norminv(p_w_hit)-norminv(p_w_fl);
%                     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
%                     D_w(4,j-4)=d_w;
%                     D_nw(4,j-4)=d_nw;
%                     Y_con4 = Y_con4+y;
%                 end
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
            Orignial=zeros(5,2);
            Switch=zeros(5,2);
            for k = 1 : 5
                Orignial(k,:) = rad(k,:)+rad(11-k,:);
                Switch (k,:)= tan(k,:)+tan(11-k,:);
            end
            ra = Orignial(:,1)./(Orignial(:,1)+Orignial(:,2));
            ta = Switch(:,1)./(Switch(:,1)+Switch(:,2));
            Y_total = [ra,ta];
            T_y = [T_y;Y_total];
        end
        
%         Orignial=zeros(5,2);
%         Switch=zeros(5,2);
%         for k = 1 : 5
%             Orignial(k,:) = rad(k,:)+rad(11-k,:);
%             Switch (k,:)= tan(k,:)+tan(11-k,:);
%         end
%         ra = Orignial(:,1)./(Orignial(:,1)+Orignial(:,2));
%         ta = Switch(:,1)./(Switch(:,1)+Switch(:,2));
%         Y_total = [ra,ta];
%         T_y = [T_y;Y_total];
    end
%     Orignial=zeros(5,2);
%     Switch=zeros(5,2);
%     for k = 1 : 5
%         Orignial(k,:) = rad(k,:)+rad(11-k,:);
%         Switch (k,:)= tan(k,:)+tan(11-k,:);
%     end
%     ra = Orignial(:,1)./(Orignial(:,1)+Orignial(:,2));
%     ta = Switch(:,1)./(Switch(:,1)+Switch(:,2));
%     Y_total = [ra,ta];
%     y_con1 = (Y_rad_con1([1,2,3,4,5],1)+Y_rad_con1([10,9,8,7,6],1))/640;
%     y_con2 = (Y_rad_con2([1,2,3,4,5],1)+Y_rad_con2([10,9,8,7,6],1))/640;
%     y_con3 = (Y_con3([1,2,3,4,5],1)+Y_con3([10,9,8,7,6],1))/640;
%     y_con4 = (Y_con4([1,2,3,4,5],1)+Y_con4([10,9,8,7,6],1))/640;   
%     yy = [y_con1;y_con2;y_con3;y_con4];
%     y_stat=[y_stat;yy];
%     y_output = [y_output;yy'];
    %% for all subject
%     y_con1 = (Y_con1([1,2,3,4,5],1)+Y_con1([10,9,8,7,6],1))/640;
%     y_con2 = (Y_con2([1,2,3,4,5],1)+Y_con2([10,9,8,7,6],1))/640;
%     y_con3 = (Y_con3([1,2,3,4,5],1)+Y_con3([10,9,8,7,6],1))/640;
%     y_con4 = (Y_con4([1,2,3,4,5],1)+Y_con4([10,9,8,7,6],1))/640;
%     y_cond1(:,kk) = (Y_con1([1,2,3,4,5],1)+Y_con1([10,9,8,7,6],1))/640;
%     y_cond2(:,kk) = (Y_con2([1,2,3,4,5],1)+Y_con2([10,9,8,7,6],1))/640;
%     y_cond3(:,kk) = (Y_con3([1,2,3,4,5],1)+Y_con3([10,9,8,7,6],1))/640;
%     y_cond4(:,kk) = (Y_con4([1,2,3,4,5],1)+Y_con4([10,9,8,7,6],1))/640;

%     Y_t=Y_con1+Y_con2+Y_con3+Y_con4;
% 
%     Y_t_plus = zeros(5,2);
%     for ll = 1 : 5 
%         Y_t_plus(ll,:) = Y_t(ll,:)+ Y_t(11-ll,:);
%     end
%     figure
%     X = categorical({'8','4','2','1','0.5'});
%     X = reordercats(X,{'8','4','2','1','0.5'});
%     %y_ax= round(Y(:,1)/(Y(1,1)+Y(1,2)) * 10000)/100 ;
%     bar(X,Y_t_plus(:,1))
% 
%     title([subject{kk},', % of trials with MS during stimulus interval ([0,50])'])
%     xlabel('tilt angle')
%     ylabel('# of trials')
%     figpath='C:\Users\86186\Desktop\fig';
%     figname=[subject{kk},'_100.png'];
%     saveas(gca, sprintf('%s\%s',figpath,figname))
end
%% run
sub_stat = [ones(40,1);ones(40,1)*2;ones(40,1)*3;ones(40,1)*4;ones(40,1)*5;ones(40,1)*6;ones(40,1)*7;ones(40,1)*8];
fa1 = [8;4;2;1;0.5];
fa1_stat = [];
for aa = 1 :64
   fa1_stat = [fa1_stat;fa1];
end
% fa2=[ones(5,1);ones(5,1)*2;ones(5,1)*3;ones(5,1)*4];
% fa2_stat = [];
% for bb = 1 : 6 
%     fa2_stat = [fa2_stat;fa2];
% end
% fact_name = {'trial difficulty','stimulus condition'};
% stats = rm_anova2(y_stat,sub_stat,fa1_stat,fa2_stat,fact_name);
amax=max(max(T_y))-min(min(T_y));
T_y = (T_y-min(min(T_y)))./amax;
mat = [T_y,sub_stat,fa1_stat];
T=array2table(mat,...
    'VariableNames',{'Radial Rate','Tan Rate','subNO','difference'});
%% 
% m1 = mat(1:20,:);
% T1 = array2table(m1,...
%     'VariableNames',{'Rate','subNO','difference','Condition'});
% 
% m2 = mat(21:40,:);
% T2 = array2table(m2,...
%     'VariableNames',{'Rate','subNO','difference','Condition'});
% 
% m3 = mat(41:60,:);
% T3 = array2table(m3,...
%     'VariableNames',{'Rate','subNO','difference','Condition'});
% 
% m4 = mat(61:80,:);
% T4 = array2table(m4,...
%     'VariableNames',{'Rate','subNO','difference','Condition'});
% 
% m5 = mat(81:100,:);
% T5 = array2table(m5,...
%     'VariableNames',{'Rate','subNO','difference','Condition'});
% 
% m6 = mat(101:120,:);
% T6 = array2table(m6,...
%     'VariableNames',{'Rate','subNO','difference','Condition'});
% %lme = fitlme(T,'Rate~1 + difference+(1|difference)+(subNO|difference)');
% % sub = categories('subNO');
% % writematrix(y_output,'A.csv'); 
% 
mdl_sum_r = fitlm(mat(:,4),mat(:,1));
figure
plot(mdl_sum_r)
xlabel('Difficult')
ylabel('MS rate')
title('RADIAL (all)')
figpath='C:\Users\86186\Desktop\fig';
name = [figpath,'/','_All_rad_reg.png'];
saveas(gca,name)

mdl_sum_t = fitlm(mat(:,4),mat(:,2));
figure
plot(mdl_sum_t)
xlabel('Difficult')
ylabel('MS rate')
title('TANGENTIAL (all)')
figpath='C:\Users\86186\Desktop\fig';
name = [figpath,'/','_All_tan_reg.png'];
saveas(gca,name)

rad_val_x1 = zeros(6,1);
tan_val_x1 = zeros(6,1);
rad_val_int = zeros(6,1);
tan_val_int = zeros(6,1);
for i = 1 : 8 
    m= mat(mat(:,3)==i,:);
%     t= tan(tan(:,2)==i,:);
    rad_all = fitlm(m(:,4),m(:,1));
    tan_all = fitlm(m(:,4),m(:,2));
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
    titlename = ['Radial Sub',num2str(i)];
    title(titlename)
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/','rad',num2str(i),'_reg.png'];
    saveas(gca,name)

    figure
    %subplot(2,1,2)
    plot(tan_all)
    xlabel('Difficult')
    ylabel('MS rate')
    titlename = ['Tangential Sub',num2str(i)];
    title(titlename)
    figpath='C:\Users\86186\Desktop\fig';
    name = [figpath,'/','tan',num2str(i),'_reg.png'];
    saveas(gca,name)
end

[h_radv0,p_radv0] = ttest(rad_val_x1);
[h_tanv0,p_tanv0] = ttest(tan_val_x1);
[h_tanvrad,p_tanvrad] = ttest(tan_val_x1,rad_val_x1);

output= [rad_val_x1,tan_val_x1];
% writematrix(output,'RvT.csv'); 
% mdl_s1 = fitlm(m1(:,3),m1(:,1));
% figure
% plot(mdl_s1)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Sub01')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','_sub1_reg.png'];
% saveas(gca,name)
% 
% mdl_s2 = fitlm(m2(:,3),m2(:,1));
% figure
% plot(mdl_s2)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Sub02')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','_sub2_reg.png'];
% saveas(gca,name)
% 
% mdl_s3 = fitlm(m3(:,3),m3(:,1));
% figure
% plot(mdl_s3)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Sub03')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','_sub3_reg.png'];
% saveas(gca,name)
% 
% mdl_s4 = fitlm(m4(:,3),m4(:,1));
% figure
% plot(mdl_s4)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Sub04')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','_sub4_reg.png'];
% saveas(gca,name)
% 
% mdl_s5 = fitlm(m5(:,3),m5(:,1));
% figure
% plot(mdl_s5)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Sub05')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','_sub5_reg.png'];
% saveas(gca,name)
% 
% mdl_s6 = fitlm(m6(:,3),m6(:,1));
% figure
% plot(mdl_s6)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Sub06')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','_sub6_reg.png'];
% saveas(gca,name)
% 
% %% RADIAL vs TANGENTIAL
% rad = mat(mat(:,4)<3,:);
% tan = mat(mat(:,4)>2,:);
% rad_all = fitlm(rad(:,3),rad(:,1));
% tan_all = fitlm(tan(:,3),tan(:,1));
% 
% figure
% %subplot(2,1,1)
% plot(rad_all)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Radial')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','radall_reg.png'];
% saveas(gca,name)
% 
% figure
% %subplot(2,1,2)
% plot(tan_all)
% xlabel('Difficult')
% ylabel('MS rate')
% title('Tangentia')
% figpath='C:\Users\86186\Desktop\fig';
% name = [figpath,'/','tanall_reg.png'];
% saveas(gca,name)
% 
% for i = 1 : 6 
%     r= rad(rad(:,2)==i,:);
%     t= tan(tan(:,2)==i,:);
%     rad_all = fitlm(r(:,3),r(:,1));
%     tan_all = fitlm(t(:,3),t(:,1));
%     figure
%     %subplot(2,1,1)
%     plot(rad_all)
%     xlabel('Difficult')
%     ylabel('MS rate')
%     titlename = ['Radial Sub',num2str(i)];
%     title(titlename)
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/','rad',num2str(i),'_reg.png'];
%     saveas(gca,name)
% 
%     figure
%     %subplot(2,1,2)
%     plot(tan_all)
%     xlabel('Difficult')
%     ylabel('MS rate')
%     titlename = ['Tangential Sub',num2str(i)];
%     title(titlename)
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/','tan',num2str(i),'_reg.png'];
%     saveas(gca,name)
% end
% 
% %% cardinal vs oblique (Cartesian)
% rad = mat(mat(:,4)==3 | mat(:,4)==1,:);
% tan = mat(mat(:,4)==2 | mat(:,4)==4,:);
% rad_all = fitlm(rad(:,3),rad(:,1));
% tan_all = fitlm(tan(:,3),tan(:,1));
% 
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
% 
% 
% for i = 1 : 6 
%     r= rad(rad(:,2)==i,:);
%     t= tan(tan(:,2)==i,:);
%     rad_all = fitlm(r(:,3),r(:,1));
%     tan_all = fitlm(t(:,3),t(:,1));
%     figure
%     %subplot(2,1,1)
%     plot(rad_all)
%     xlabel('Difficult')
%     ylabel('MS rate')
%     titlename = ['Cardinal Sub',num2str(i)];
%     title(titlename)
%     figpath='C:\Users\86186\Desktop\fig';
%     name = [figpath,'/','Cardinal_Cartesian',num2str(i),'_reg.png'];
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
%     name = [figpath,'/','Oblique_Cartesian',num2str(i),'_reg.png'];
%     saveas(gca,name)
% end