subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};    
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
        for j = 1 : 8

            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{1}, ...
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
            if i == 2 
                if j < 5
                    %acc
                    ACC_W(1,j)=mean(acc_w);
                    ACC_NW(1,j)=mean(acc_nw);
                    w_num_con1 = w_num_con1+size(acc_w,1);
                    nw_num_con1= nw_num_con1+size(acc_nw,1);
                    %rec
                    REC_W(1,j)=mean(rec_w(rec_w<4));
                    REC_NW(1,j)=mean(rec_nw(rec_nw<4));
                    %d_prime
                    p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
                    p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
                    p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
                    p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
                    d_w=norminv(p_w_hit)-norminv(p_w_fl);
                    d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
                    D_w(1,j)=d_w;
                    D_nw(1,j)=d_nw;
                    
                else
                    ACC_W(2,j-4)=mean(acc_w);
                    ACC_NW(2,j-4)=mean(acc_nw);
                    w_num_con2 = w_num_con2+size(acc_w,1);
                    nw_num_con2= nw_num_con2+size(acc_nw,1);
                    %rec
                    REC_W(2,j-4)=mean(rec_w(rec_w<4));
                    REC_NW(2,j-4)=mean(rec_nw(rec_nw<4));
                    %d_prime
                    p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
                    p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
                    p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
                    p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
                    d_w=norminv(p_w_hit)-norminv(p_w_fl);
                    d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
                    D_w(2,j-4)=d_w;
                    D_nw(2,j-4)=d_nw;
                    
                end
            else
                if j < 5
                    ACC_W(3,j)=mean(acc_w);
                    ACC_NW(3,j)=mean(acc_nw);
                    w_num_con3 = w_num_con3+size(acc_w,1);
                    nw_num_con3= nw_num_con3+size(acc_nw,1);
                    %rec
                    REC_W(3,j)=mean(rec_w(rec_w<4));
                    REC_NW(3,j)=mean(rec_nw(rec_nw<4));
                    %d_prime
                    p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
                    p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
                    p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
                    p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
                    d_w=norminv(p_w_hit)-norminv(p_w_fl);
                    d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
                    D_w(3,j)=d_w;
                    D_nw(3,j)=d_nw;
                    
                else
                    ACC_W(4,j-4)=mean(acc_w);
                    ACC_NW(4,j-4)=mean(acc_nw);
                    w_num_con4 = w_num_con4+size(acc_w,1);
                    nw_num_con4= nw_num_con4+size(acc_nw,1);
                    %rec
                    REC_W(4,j-4)=mean(rec_w(rec_w<4));
                    REC_NW(4,j-4)=mean(rec_nw(rec_nw<4));
                    %d_prime
                    p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
                    p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
                    p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
                    p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
                    d_w=norminv(p_w_hit)-norminv(p_w_fl);
                    d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
                    D_w(4,j-4)=d_w;
                    D_nw(4,j-4)=d_nw;
                end
            end
        end
    end

    
%data: 4组，每组5个数据
% y=[0.34 0.32 0.34 0.47 0.51; 
%    0.44 0.50 0.53 0.56 0.61;
%    0.62 0.52 0.56 0.92 0.94;
%    0.46 0.44 0.47 0.65 0.68;
% ];


acc_w1=mean(ACC_W,2);
acc_nw1=mean(ACC_NW,2);
ACC_TOTAL=[acc_w1,acc_nw1];
% acc_w1_sem=std(ACC_W,2)/2;
% acc_nw1_sem=std(ACC_NW,2)/2;
% ACC_TOTAL_SEM=[acc_w1_sem,acc_nw1_sem];


x=[1 2 3 4]; %对应4组
t1=ACC_TOTAL(:,1); 
t2=ACC_TOTAL(:,2); 
% t3=y(:,3);
% t4=y(:,4);
% t5=y(:,5);
t=[t1;t2];
% er=errorbar(X,ACC_TOTAL,ACC_TOTAL_SEM);
% er.LineStyle= 'none';  
%plot bar
c=bar(x,ACC_TOTAL);
set(c(1,1),'FaceColor','m','BarWidth',0.9);% c(1,1) 就是设置第一组的第一个数据柱
set(c(1,2),'FaceColor','Y','BarWidth',0.9);% c(1,2) 设置第二个
% set(c(1,3),'FaceColor','R','BarWidth',0.9);%‘FaceColor’ 设置柱的颜色
% set(c(1,4),'FaceColor','G','BarWidth',0.9);%‘BarWidth’ 设置柱的宽度
% set(c(1,5),'FaceColor','B','BarWidth',0.9);% 对于单个柱，设置c(1),c(2),c(3),...即可


num_w = [w_num_con1,w_num_con2,w_num_con3,w_num_con4];
num_nw=[nw_num_con1,nw_num_con2,nw_num_con3,nw_num_con4];
% 在每个柱顶部添加相应的数字
for i = 1:4    
    text(x(i)-0.15,t1(i),num2str(num_w(i)),...   
        'HorizontalAlignment','center',...    
        'VerticalAlignment','bottom','FontSize',10)  %0.3可以调整相邻柱之间的间隔，手动调节 
    text(x(i)+0.15,t2(i),num2str(num_nw(i)),...    
        'HorizontalAlignment','center',...    
        'VerticalAlignment','bottom','FontSize',10)   %‘FontSize’设置文字大小
%     text(x(i),t3(i),num2str(t3(i)),...   
%         'HorizontalAlignment','center',...   
%         'VerticalAlignment','bottom','FontSize',15)
%     text(x(i)+0.15,t4(i),num2str(t4(i)),...   
%         'HorizontalAlignment','center',...   
%         'VerticalAlignment','bottom','FontSize',15)
%     text(x(i)+0.3,t5(i),num2str(t5(i)),...   
%         'HorizontalAlignment','center',...   
%         'VerticalAlignment','bottom','FontSize',15)
end
 
%设置 x y轴刻度标签
set(gca,'XTickLabel',{'Cardinal directions ON meridians','Oblique directions OFF meridians','Cardinal directions OFF meridians','Oblique directions ON meridians'},'FontSize',7);
set(gca,'YTicklabel',{'0','0.2','0.4','0.6','0.8','1'});
 
%设置 y轴刻度
set(gca,'YTick',0:0.2:1);
% set(gca,'legend',{'a','b'}); 
% set(gca,'title',{'accuracy'}); 
%给 x y轴加名字
% xlabel('condition','FontSize',15);
ylabel('% correct','FontSize',10);
% legend('a','b');
title('accuracy');
legend('aaa','bbb');
saveas(gca, sprintf('%s/a.png',MATpath)) 
%不同算法图例
%legend('A','B','C','D','Ours','location','northwest');%'location' 控制图例放置位置
 
% %将结果图像保存为PDF格式
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% filename='result';%只需改动名字
% print(gcf,filename,'-dpdf','-r0')
% close(gcf);
%  
% %若保存的PDF有较多空白边缘，则修改纸的大小如下
% fig=gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3)-5 fig_pos(4)-4];%PaperSize: WxH 可以调节纸大小，进而修改空白边
% filename='result';
% print(gcf,filename,'-dpdf','-r0');