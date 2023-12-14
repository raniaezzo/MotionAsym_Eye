clc; clear all;

tab_total = [] ;
MS_total = [] ; 
sample=[];
sub={'S01','S02','S03','S04','S05','S06','S07','S08'};
subject = 'S01'; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
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
figpath='C:\Users\86186\Desktop\fig\new';

color={[215,25,28],[253,174,97],[255,255,191],[171,217,233],[44,123,182]};	

for k = 1 : 8
    Y=zeros(10,2);
    y_total=0;
    for i = 1 : 2
        if i ==1 & k > 6
            continue
        end
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
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction_name{j}) );
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
    bs = bar(X,y_ax);
    bs.FaceColor = 'flat';
    for sd = 1 : 5 
        bs.CData(sd,:) = color{6-sd}/255;
        bs.CData((11-sd),:) = color{6-sd}/255;
    end
    
    
    title_na = ['Example Subject ',sub{k}];
    title(title_na)
    xlabel('tilt from standard (degrees)')
    ylabel('% trials with microsaccade during stimulus')
    ylim([0,ceil(max(y_ax))+2])
    for s = 1 : 10
        str1 = [num2str(y_ax(s)),'%'];
        text(s,y_ax(s),str1,...   
               'HorizontalAlignment','center',...    
               'VerticalAlignment','bottom')
    end  
    name = [figpath,'/',sub{k},'-tilt2.png'];
%     title(title_name)
    saveas(gca,name)
    
end
