% % count different session
% clc; clear all;
% 
% tab_total = [] ;
% MS_total = [] ; 
% sample=[];
% subject = 'S08'; condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'}; 
% direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
% ACC_W=[];
% ACC_NW=[];
% REC_W=[];
% REC_NW=[];
% Y=zeros(10,2);
% Tab_w_cl=[];
% Tab_nw_cl=[];
% Tab_w_ncl=[];
% Tab_nw_ncl=[];
% 
% 
% for i = 2 : 2
%     for j = 1 : 8
%         
%         main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject, ...
%             'RawData', condition{i}, 'Block1');
%         cd(fullfile(main_folder, 'eyedata'));
%         edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
%         edf_path = fullfile(main_folder,'eyedata',edf_name);
%         msg_filepath=replace(edf_path,'edf','msg');
%         samplingRateData=findSamplingRate(msg_filepath);
%         sample=[sample,samplingRateData];
%         MATpath = fullfile(main_folder, 'eyedata','MATs');
% %         cd(fullfile(main_folder, 'eyedata','MATs'));
% %         ms_name = dir(sprintf('%s.mat', direction{j})).name;
%         ms_path= fullfile(MATpath,sprintf('%s.mat', direction{j}) );
%         load(ms_path);
% %         MS_TEMP(:,10)=MS_TEMP(:,10)+ ((i-1)*8 + (j-1))*800;   
% %         MS_total=[MS_total;MS_TEMP];
% %         MATpath = fullfile(main_folder, 'eyedata','MATs');
%         tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
%         load(tab_path)
% %         tab(:,1)=tab(:,1)+ ((i-1)*8 + (j-1))*800;  
% %         tab_total=[tab_total;tab];
%         [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,20);
%         ACC_W=[ACC_W;acc_w];
%         ACC_NW=[ACC_NW;acc_nw];
%         REC_W=[REC_W;rec_w];
%         REC_NW=[REC_NW;rec_nw];
%         Y=Y+y;
%         Tab_w_cl=[Tab_w_cl;tab_w_cl];
%         Tab_nw_cl=[Tab_nw_cl;tab_nw_cl];
%         Tab_w_ncl=[Tab_w_ncl;tab_w_ncl];
%         Tab_nw_ncl=[Tab_nw_ncl;tab_nw_ncl];
%     end
%     
%     
% end
% p_w_hit=sum(Tab_w_cl(:,15))/size(Tab_w_cl,1);
% p_nw_hit=sum(Tab_nw_cl(:,15))/size(Tab_nw_cl,1);
% p_w_fl=(size(Tab_w_ncl,1)-sum(Tab_w_ncl(:,14)))/size(Tab_w_ncl,1);
% p_nw_fl=(size(Tab_nw_ncl,1)-sum(Tab_nw_ncl(:,14)))/size(Tab_nw_ncl,1);
% 
% d_w=norminv(p_w_hit)-norminv(p_w_fl);
% d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
% 
% x=1:10;
% 
% figure
% X = categorical({'with MS in sti','notwith MS in sti'});
% X = reordercats(X,{'with MS in sti','notwith MS in sti'});
% bar(X,[mean(ACC_W),mean(ACC_NW)])
% 
% title('Accuracy')
% saveas(gca, sprintf('%s/ACC.png',MATpath))
% 
% 
% figure
% X = categorical({'with MS in sti','notwith MS in sti'});
% X = reordercats(X,{'with MS in sti','notwith MS in sti'});
% bar(X,[mean(REC_W(REC_W<4)),mean(REC_NW(REC_NW<4))])
% hold on 
% er=errorbar(X,[mean(REC_W(REC_W<4)),mean(REC_NW(REC_NW<4))],[std(REC_W(REC_W<4)),std(REC_NW(REC_NW<4))]);
% er.LineStyle= 'none';  
% hold off 
% 
% title('reaction time')
% saveas(gca, sprintf('%s/REC.png',MATpath))
% 
% figure
% X = categorical({'with MS in sti','notwith MS in sti'});
% X = reordercats(X,{'with MS in sti','notwith MS in sti'});
% bar(X,[d_w,d_nw])
% 
% title('d-prime')
% saveas(gca, sprintf('%s/dprime.png',MATpath))
% 
% 
% figure
% X = categorical({'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
% X = reordercats(X,{'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
% bar(X,Y,'stacked')
% 
% title('stacked')
% saveas(gca, sprintf('%s/deg_stacked.png',MATpath))

    


clc; clear all;

tab_total = [] ;
MS_total = [] ; 
sample=[];
subject = 'S01'; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
ACC_W=[];
ACC_NW=[];
REC_W=[];
REC_NW=[];
Y=zeros(10,2);
Tab_w_cl=[];
Tab_nw_cl=[];
Tab_w_ncl=[];
Tab_nw_ncl=[];


for i = 1 : 2
    for j = 1 : 8
        
        main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject, ...
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
        [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,20);
        ACC_W=[ACC_W;acc_w];
        ACC_NW=[ACC_NW;acc_nw];
        REC_W=[REC_W;rec_w];
        REC_NW=[REC_NW;rec_nw];
        Y=Y+y;
        Tab_w_cl=[Tab_w_cl;tab_w_cl];
        Tab_nw_cl=[Tab_nw_cl;tab_nw_cl];
        Tab_w_ncl=[Tab_w_ncl;tab_w_ncl];
        Tab_nw_ncl=[Tab_nw_ncl;tab_nw_ncl];
    end
    
    
end
p_w_hit=sum(Tab_w_cl(:,15))/size(Tab_w_cl,1);
p_nw_hit=sum(Tab_nw_cl(:,15))/size(Tab_nw_cl,1);
p_w_fl=(size(Tab_w_ncl,1)-sum(Tab_w_ncl(:,14)))/size(Tab_w_ncl,1);
p_nw_fl=(size(Tab_nw_ncl,1)-sum(Tab_nw_ncl(:,14)))/size(Tab_nw_ncl,1);

d_w=norminv(p_w_hit)-norminv(p_w_fl);
d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);

x=1:10;

figure
X = categorical({'with MS in sti','notwith MS in sti'});
X = reordercats(X,{'with MS in sti','notwith MS in sti'});
bar(X,[mean(ACC_W),mean(ACC_NW)])

title('Accuracy')
ylabel('% correct')
saveas(gca, sprintf('%s/ACC.png',MATpath))


figure
X = categorical({'with MS in sti','notwith MS in sti'});
X = reordercats(X,{'with MS in sti','notwith MS in sti'});
bar(X,[mean(REC_W(REC_W<4)),mean(REC_NW(REC_NW<4))])
hold on 
er=errorbar(X,[mean(REC_W(REC_W<4)),mean(REC_NW(REC_NW<4))],[std(REC_W(REC_W<4)),std(REC_NW(REC_NW<4))]);
er.LineStyle= 'none';  
hold off 

title('reaction time')
ylabel('sec(s)')
saveas(gca, sprintf('%s/REC.png',MATpath))

figure
X = categorical({'with MS in sti','notwith MS in sti'});
X = reordercats(X,{'with MS in sti','notwith MS in sti'});
bar(X,[d_w,d_nw])

title('d-prime')
saveas(gca, sprintf('%s/dprime.png',MATpath))


figure
X = categorical({'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
X = reordercats(X,{'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
bar(X,Y,'stacked')

title('With vs. Without MS')
xlabel('tilt angle')
ylabel('# of trials')
saveas(gca, sprintf('%s/deg_stacked.png',MATpath))