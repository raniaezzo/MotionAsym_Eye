%% this is for the 64 conditions
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'};
conditions = {'Switched', 'original'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
location = {'VU','VL','HL','HR','LL','LR','UL','UR'};
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
acc_with_sub_car =zeros(5,6);
acc_notwith_sub_car =zeros(5,6);
acc_with_sub_obl =zeros(5,6);
acc_notwith_sub_obl =zeros(5,6);
for kk = 6 : 6
    tab_nw_total_car = [];
    tab_w_total_car = [];
    tab_nw_total_obl = [];
    tab_w_total_obl = [];
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
    acc_with_car =zeros(5,1);
    acc_notwith_car = zeros(5,1);
    acc_with_obl =zeros(5,1);
    acc_notwith_obl = zeros(5,1);
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
            num_w = zeros(1,5);
            num_nw = zeros(1,5);
            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,padding_time);
            tab_nw = [tab_nw_cl;tab_nw_ncl];
            tab_w = [tab_w_cl;tab_w_ncl];
            if tab_nw(1,9) < 5
                
                    for jj = 1 : 4
                        tab_nw = [tab_nw_cl;tab_nw_ncl];
                        tab_w = [tab_w_cl;tab_w_ncl];
                        tab_nw_r = tab_nw(tab_nw(:,9)==jj,:);
                        tab_w_r = tab_w(tab_w(:,9)==jj,:);
                        for ti = 1:5
                            a_nw_car = mean(tab_nw_r(tab_nw_r(:,11)==tilt(ti),15));
                            num_nw(ti) = size(tab_nw_r(tab_nw_r(:,11)==tilt(ti)),1);
                            a_w_car = mean(tab_w_r(tab_w_r(:,11)==tilt(ti),15));
                            num_w(ti) = size(tab_w_r(tab_w_r(:,11)==tilt(ti)),1);
                            acc_with_car(ti)=a_w_car;
                            acc_notwith_car(ti) = a_nw_car;
                        end
                        
                        num_wanw = [num_w;num_nw];            
                        a_wnw_car = [acc_with_car,acc_notwith_car]';
            %             a_wnw_obl = [acc_with_obl,acc_notwith_obl]';
                        x = 1:5;
                        figure
                        c= bar(x,a_wnw_car);
                        set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','with MS');
                        set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','without MS');
                        for ii = 1:5    
                            text(x(ii)-0.15,a_wnw_car(1,ii),num2str(round(num_wanw(1,ii),2)),...   
                                'HorizontalAlignment','center',...    
                                'VerticalAlignment','bottom','FontSize',6) 
                            text(x(ii)+0.15,a_wnw_car(2,ii),num2str(round(num_wanw(2,ii),2)),...    
                                'HorizontalAlignment','center',...    
                                'VerticalAlignment','bottom','FontSize',6)  
                        end
                        ylim([0.3,1.1])
                        legend(c(1,1:2),'With MS','Not with MS')
                        xticklabels({'8°','4°','2°','1°','0.5°'})
                        ylabel('Accuracy')
                        titlename = [subject{kk},'_ ',direction{j},'_ ',location{jj},' Accuracy'];
                        title(titlename,'FontSize',12);
                        figpath='C:\Users\86186\Desktop\fig';
                        name = [figpath,'/',subject{kk},'_ ',direction{j},'_ ',location{jj},'_acc.png'];
                        saveas(gca,name)               
                    end
            else
                for jj = 5:8
                        tab_nw = [tab_nw_cl;tab_nw_ncl];
                        tab_w = [tab_w_cl;tab_w_ncl];
                        tab_nw_r = tab_nw(tab_nw(:,9)==jj,:);
                        tab_w_r = tab_w(tab_w(:,9)==jj,:);
                        for ti = 1:5
                            a_nw_car = mean(tab_nw_r(tab_nw_r(:,11)==tilt(ti),15));
                            num_nw(ti) = size(tab_nw_r(tab_nw_r(:,11)==tilt(ti)),1);
                            a_w_car = mean(tab_w_r(tab_w_r(:,11)==tilt(ti),15));
                            num_w(ti) = size(tab_w_r(tab_w_r(:,11)==tilt(ti)),1);
                            acc_with_car(ti)=a_w_car;
                            acc_notwith_car(ti) = a_nw_car;
                        end
                        
                        num_wanw = [num_w;num_nw];            
                        a_wnw_car = [acc_with_car,acc_notwith_car]';
            %             a_wnw_obl = [acc_with_obl,acc_notwith_obl]';
                        x = 1:5;
                        figure
                        c= bar(x,a_wnw_car);
                        set(c(1,1),'FaceColor','m','BarWidth',0.9,'DisplayName','with MS');
                        set(c(1,2),'FaceColor','Y','BarWidth',0.9,'DisplayName','without MS');
                        for ii = 1:5    
                            text(x(ii)-0.15,a_wnw_car(1,ii),num2str(round(num_wanw(1,ii),2)),...   
                                'HorizontalAlignment','center',...    
                                'VerticalAlignment','bottom','FontSize',6) 
                            text(x(ii)+0.15,a_wnw_car(2,ii),num2str(round(num_wanw(2,ii),2)),...    
                                'HorizontalAlignment','center',...    
                                'VerticalAlignment','bottom','FontSize',6)  
                        end
                        ylim([0.3,1.1])
                        legend(c(1,1:2),'With MS','Not with MS')
                        xticklabels({'8°','4°','2°','1°','0.5°'})
                        ylabel('Accuracy')
                        titlename = [subject{kk},'_ ',direction{j},'_ ',location{jj},' Accuracy'];
                        title(titlename,'FontSize',12);
                        figpath='C:\Users\86186\Desktop\fig';
                        name = [figpath,'/',subject{kk},'_ ',direction{j},'_ ',location{jj},'_acc.png'];
                        saveas(gca,name)     
                end
            end
                
          

        end

    
    end




end