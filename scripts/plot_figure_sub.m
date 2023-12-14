clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};

for kk = 1 : 6 
    tab_total = [] ;
    MS_total = [] ; 
    sample=[];
    ACC_W=[];
    ACC_NW=[];
    REC_W=[];
    REC_NW=[];
%     ACC_W=zeros(4,4);
%     ACC_NW=zeros(4,4);
%     REC_W=zeros(4,4);
%     REC_NW=zeros(4,4);
    Y=zeros(10,2);
    Tab_w_cl=[];
    Tab_nw_cl=[];
    Tab_w_ncl=[];
    Tab_nw_ncl=[];
    D_w=[];
    D_nw=[];
    w_num=0;
    nw_num=0;
%     D_w=zeros(4,4);
%     D_nw=zeros(4,4);
   %calculate numbers of different condition
%     w_num_con1 = 0;
%     nw_num_con1= 0;
%     w_num_con2 = 0;
%     nw_num_con2= 0;
%     w_num_con3 = 0;
%     nw_num_con3= 0;
%     w_num_con4 = 0;
%     nw_num_con4= 0;
    
    
    padding_time=[0,50];
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
            Y=Y+y;
%             if i == 2 
%                 if j < 5
%                     %acc
%                     ACC_W(1,j)=mean(acc_w);
%                     ACC_NW(1,j)=mean(acc_nw);
%                     w_num_con1 = w_num_con1+size(acc_w,1);
%                     nw_num_con1= nw_num_con1+size(acc_nw,1);
%                     %rec
%                     REC_W(1,j)=mean(rec_w(rec_w<4));
%                     REC_NW(1,j)=mean(rec_nw(rec_nw<4));
%                     %d_prime
%                     p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
%                     p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
%                     p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
%                     p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
%                     d_w=norminv(p_w_hit)-norminv(p_w_fl);
%                     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
%                     D_w(1,j)=d_w;
%                     D_nw(1,j)=d_nw;
%                     
%                 else
%                     ACC_W(2,j-4)=mean(acc_w);
%                     ACC_NW(2,j-4)=mean(acc_nw);
%                     w_num_con2 = w_num_con2+size(acc_w,1);
%                     nw_num_con2= nw_num_con2+size(acc_nw,1);
%                     %rec
%                     REC_W(2,j-4)=mean(rec_w(rec_w<4));
%                     REC_NW(2,j-4)=mean(rec_nw(rec_nw<4));
%                     %d_prime
%                     p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
%                     p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
%                     p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
%                     p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);
%                     d_w=norminv(p_w_hit)-norminv(p_w_fl);
%                     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);
%                     D_w(2,j-4)=d_w;
%                     D_nw(2,j-4)=d_nw;
%                     
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
%                     
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
%                 end
%             end
          
            ACC_W=[ACC_W;mean(acc_w)];
            w_num = w_num+size(acc_w,1);
            ACC_NW=[ACC_NW;mean(acc_nw)];
            nw_num = nw_num+size(acc_nw,1);

            REC_W=[REC_W;mean(rec_w(rec_w<4))];
            REC_NW=[REC_NW;mean(rec_nw(rec_nw<4))];
%             Y=Y+y;
            p_w_hit=sum(tab_w_cl(:,15))/size(tab_w_cl,1);
            p_nw_hit=sum(tab_nw_cl(:,15))/size(tab_nw_cl,1);
            p_w_fl=(size(tab_w_ncl,1)-sum(tab_w_ncl(:,14)))/size(tab_w_ncl,1);
            p_nw_fl=(size(tab_nw_ncl,1)-sum(tab_nw_ncl(:,14)))/size(tab_nw_ncl,1);

            d_w=norminv(p_w_hit)-norminv(p_w_fl);
            d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);

            D_w=[D_w;d_w];
            D_nw=[D_nw;d_nw];
%             Tab_w_cl=[Tab_w_cl;tab_w_cl];
%             Tab_nw_cl=[Tab_nw_cl;tab_nw_cl];
%             Tab_w_ncl=[Tab_w_ncl;tab_w_ncl];
%             Tab_nw_ncl=[Tab_nw_ncl;tab_nw_ncl];
        end


    end
%     total_num = [w_num,nw_num];
%     p_w_hit=sum(Tab_w_cl(:,15))/size(Tab_w_cl,1);
%     p_nw_hit=sum(Tab_nw_cl(:,15))/size(Tab_nw_cl,1);
%     p_w_fl=(size(Tab_w_ncl,1)-sum(Tab_w_ncl(:,14)))/size(Tab_w_ncl,1);
%     p_nw_fl=(size(Tab_nw_ncl,1)-sum(Tab_nw_ncl(:,14)))/size(Tab_nw_ncl,1);
%     
%     d_w=norminv(p_w_hit)-norminv(p_w_fl);
%     d_nw=norminv(p_nw_hit)-norminv(p_nw_fl);

    x=1:10;

%     figure
%     X = categorical({'with MS in sti','notwith MS in sti'});
%     X = reordercats(X,{'with MS in sti','notwith MS in sti'});
%     bar(X,[mean(ACC_W),mean(ACC_NW)])
%     hold on 
%     st_acc=std(ACC_W);
%     sem_acc=st_acc/sqrt(16);
%     st_nacc=std(ACC_NW);
%     sem_nacc=st_nacc/sqrt(16);
% 
%     y_axis=[mean(ACC_W),mean(ACC_NW)];
% 
% 
%     er=errorbar(X,[mean(ACC_W),mean(ACC_NW)],[sem_acc,sem_nacc]);
%     er.LineStyle= 'none';  
%     hold off 
%     title('Accuracy')
%     ylabel('% correct')
%     for i = 1 : 2
%         str = ['total number :',num2str(total_num(i))];
%         text(i,y_axis(i)/2,str,...   
%                 'HorizontalAlignment','center',...    
%                 'VerticalAlignment','bottom','FontSize',15)
%     end
%     saveas(gca, sprintf('%s/ACC2.png',MATpath))
% 
% 
%     figure
%     X = categorical({'with MS in sti','notwith MS in sti'});
%     X = reordercats(X,{'with MS in sti','notwith MS in sti'});
%     bar(X,[mean(REC_W),mean(REC_NW)])
%     hold on 
%     st_rec=std(REC_W);
%     sem_rec=st_rec/sqrt(16);
%     st_nrec=std(REC_NW);
%     sem_nrec=st_nrec/sqrt(16);
% 
%     er=errorbar(X,[mean(REC_W),mean(REC_NW)],[sem_rec,sem_nrec]);
%     er.LineStyle= 'none';  
%     hold off 
% 
%     title('reaction time')
%     ylabel('sec(s)')
% 
%     for i = 1 : 2
%         str = ['total number :',num2str(total_num(i))];
%         text(i,y_axis(i)/2,str,...   
%                 'HorizontalAlignment','center',...    
%                 'VerticalAlignment','bottom','FontSize',15)
%     end
% 
%     saveas(gca, sprintf('%s/REC2.png',MATpath))

%     figure
%     X = categorical({'with MS in sti','notwith MS in sti'});
%     X = reordercats(X,{'with MS in sti','notwith MS in sti'});
%     bar(X,[mean(D_w(D_w<10)),mean(D_nw(D_nw<10))])
%     hold on 
%     st_dp=std(D_w(D_w<10));
%     sem_dp=st_dp/sqrt(16);
%     st_ndp=std(D_nw(D_nw<10));
%     sem_ndp=st_ndp/sqrt(16);
% 
%     er=errorbar(X,[mean(D_w(D_w<10)),mean(D_nw(D_nw<10))],[sem_dp,sem_ndp]);
%     er.LineStyle= 'none';  
%     hold off 
% 
%     title('d-prime')
%     saveas(gca, sprintf('%s/dprime2.png',MATpath))


    figure
    X = categorical({'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
    X = reordercats(X,{'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
    y_ax= round(Y(:,1)/(Y(1,1)+Y(1,2)) * 10000)/100 ;
    bar(X,y_ax)

    title('% of trials with MS during stimulus interval')
    xlabel('tilt angle')
    ylabel('% of trials')

    for k = 1 : 10
        str1 = [num2str(y_ax(k)),'%'];
        text(k,y_ax(k),str1,...   
                'HorizontalAlignment','center',...    
                'VerticalAlignment','bottom')
    end

    saveas(gca, sprintf('%s/deg_stacked.png',MATpath))
end