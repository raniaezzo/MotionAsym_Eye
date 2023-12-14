%Cardinal vs. Oblique
%% this script is changed by plot_figure_diffcon2.m and this is for plot MS count figure
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
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

s1 = [];
s2 = [];
figpath='C:\Users\86186\Desktop\fig';
for kk = 1 : 8 
%     Y_con1 = 0 ; 
%     Y_con2 = 0 ; 
%     Y_con3 = 0 ; 
%     Y_con4 = 0 ; 
%     Orig = zeros(10,2);
%     Swit = zeros(10,2);
    Card = zeros(10,2);
    Obli = zeros(10,2);
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
    for i = 1 : 2
        if kk > 6 & i == 1 
            continue
        end
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
            if j < 5
               Card = Card + y; 
            else
               Obli = Obli + y;
                
            end
          

        end

    
    end
    Orignial=zeros(5,2);
    Switch=zeros(5,2);
    for k = 1 : 5
        Orignial(k,:) = Card(k,:)+Card(11-k,:);
        Switch (k,:)= Obli(k,:)+Obli(11-k,:);
    end
    
    
%     ca = Orignial(:,1)./(Orignial(:,1)+Orignial(:,2));
%     ob = Switch(:,1)./(Switch(:,1)+Switch(:,2));
    ca = Orignial(:,1);
    ob = Switch(:,1);
%    Y_total = [ca,ob]/1280;
    s1 = [s1;sum(ca)];
    s2 = [s2;sum(ob)];
    
    

%    saveas(gca,name)


    
       
end
color = {[127, 191, 123]	,[175, 141, 195]}	;	
figpath_new = 'C:\Users\86186\Desktop\fig\new';
ss = [s1,s2];
X = categorical({'Cardinal','Oblique'});
X = reordercats(X,{'Cardinal','Oblique'});
x=[1,2];
for sdf = 1 : 8
    if sdf < 7         
        figure
        c = bar(X,ss(sdf,:)/6400);
        c.FaceColor = 'flat';
        c.CData(1,:) = color{1}/255;
        c.CData(2,:) = color{2}/255;
        title_na = ['Example Subject ',subject{sdf}];
        title(title_na)
        ylabel('% trials with microsaccade during stimulus')
        name = [figpath_new,'/',subject{sdf},'-CvO.png'];
        saveas(gca,name)
    else
        figure
        c = bar(X,ss(sdf,:)/3200);
        c.FaceColor = 'flat';
        c.CData(1,:) = color{1}/255;
        c.CData(2,:) = color{2}/255;
        title_na = ['Example Subject ',subject{sdf}];
        title(title_na)
        ylabel('% trials with microsaccade during stimulus')
        name = [figpath_new,'/',subject{sdf},'-CvO.png'];
        saveas(gca,name)
    end
        
end
    

ss_nor = zeros(8,2);
ss_nor(1:6,:) = ss(1:6,:)/6400;
ss_nor([7,8],:) = ss([7,8],:)/6400;
ss_nor_plot = mean(ss_nor,1);
figure
c = bar(X,ss_nor_plot);
c.FaceColor = 'flat';
c.CData(1,:) = color{1}/255;
c.CData(2,:) = color{2}/255;
%title_na = ['Example Subject ',subject{sdf}];
title('% of trial with microsaccade during stimulus')
ylabel('% of trials')
name = 'C:\Users\86186\Desktop\fig\new\tt.png';
saveas(gca,name)



% figure
% c = bar(X,[s1,s2],'stacked');
% %c.FaceColor = 'flat';
% color = {[51, 34, 136],[17, 119, 51],[68, 170, 153],[136, 204, 238],[221, 204, 119],[204, 102, 119],[170, 68, 153],[136, 34 85]};
% for sdf = 1 : 8 
%     c(sdf).FaceColor = 'flat';
%     c(sdf).CData = color{sdf}/255;
% end
% 
% ylabel('# trials with microsaccade during stimulus')
% title('All subject')
% name = 'C:\Users\86186\Desktop\fig\new\CvO.png';
%saveas(gca,name)