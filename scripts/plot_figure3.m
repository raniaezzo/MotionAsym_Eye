clc; clear all;

tab_total = [] ;
MS_total = [] ; 
sample=[];
sub={'S01','S02','S03','S04','S05','S06','S07','S08'};
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

acc_total_w=[];
acc_total_nw=[];
rec_total_w=[];
rec_total_nw=[];
d_total_w=[];
d_total_nw=[];


for k = 1 : 1
    y_total=0;
    for i = 2 : 2
        for j = 1 : 8

            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', sub{k}, ...
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
            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,[0,100]);

            Y=Y+y;

        end


    end

    y_total=y_total+Y;
    
    figure
    X = categorical({'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
    X = reordercats(X,{'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
    % bar(X,Y,'stacked')
    y_ax= round(y_total(:,1)/(y_total(1,1)+y_total(1,2)) * 10000)/100 ;
    bar(X,y_ax)
    title('With vs. Without MS')
    xlabel('tilt from standard (degrees)')
    ylabel('% trials with microsaccade during stimulus')
    ylim([0,ceil(max(y_ax))+1])
    for s = 1 : 10
        str1 = [num2str(y_ax(s)),'%'];
        text(s,y_ax(s),str1,...   
               'HorizontalAlignment','center',...    
               'VerticalAlignment','bottom')
    end    
    
end


% x=1:10;
% 
% figure
% X = categorical({'with MS in sti','notwith MS in sti'});
% X = reordercats(X,{'with MS in sti','notwith MS in sti'});
% bar(X,[mean(ACC_W),mean(ACC_NW)])
% hold on 
% st_acc=std(ACC_W);
% sem_acc=st_acc/sqrt(7);
% st_nacc=std(ACC_NW);
% sem_nacc=st_nacc/sqrt(7);
% 
% er=errorbar(X,[mean(ACC_W),mean(ACC_NW)],[sem_acc,sem_nacc]);
% er.LineStyle= 'none';  
% hold off 
% title('Accuracy')
% ylabel('% correct')
% saveas(gca, sprintf('%s/ACC2.png',MATpath))
% 
% 
% figure
% X = categorical({'with MS in sti','notwith MS in sti'});
% X = reordercats(X,{'with MS in sti','notwith MS in sti'});
% bar(X,[mean(REC_W),mean(REC_NW)])
% hold on 
% st_rec=std(REC_W);
% sem_rec=st_rec/sqrt(7);
% st_nrec=std(REC_NW);
% sem_nrec=st_nrec/sqrt(7);
% 
% er=errorbar(X,[mean(REC_W),mean(REC_NW)],[sem_rec,sem_nrec]);
% er.LineStyle= 'none';  
% hold off 
% 
% title('reaction time')
% ylabel('sec(s)')
% saveas(gca, sprintf('%s/REC2.png',MATpath))
% 
% figure
% X = categorical({'with MS in sti','notwith MS in sti'});
% X = reordercats(X,{'with MS in sti','notwith MS in sti'});
% bar(X,[mean(D_w(D_w<10)),mean(D_nw(D_nw<10))])
% hold on 
% st_dp=std(D_w(D_w<10));
% sem_dp=st_dp/sqrt(7);
% st_ndp=std(D_nw(D_nw<10));
% sem_ndp=st_ndp/sqrt(7);
% 
% er=errorbar(X,[mean(D_w(D_w<10)),mean(D_nw(D_nw<10))],[sem_dp,sem_ndp]);
% er.LineStyle= 'none';  
% hold off 
% 
% title('d-prime')
% saveas(gca, sprintf('%s/dprime2.png',MATpath))
% 
% 
% figure
% X = categorical({'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
% X = reordercats(X,{'-8','-4','-2','-1','-0.5','0.5','1','2','4','8'});
% % bar(X,Y,'stacked')
% y_ax= round(Y(:,1)/(Y(1,1)+Y(1,2)) * 10000)/100 ;
% bar(X,y_ax)
% title('With vs. Without MS')
% xlabel('tilt angle')
% ylabel('% of trials')
% 
% for k = 1 : 10
%     str1 = [num2str(y_ax(k)),'%'];
%     text(k,y_ax(k),str1,...   
%            'HorizontalAlignment','center',...    
%            'VerticalAlignment','bottom')
% end
% saveas(gca, sprintf('%s/deg_stacked.png',MATpath))