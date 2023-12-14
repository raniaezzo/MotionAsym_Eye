%% same location/ stimulus MS direction
clear all; clc;

directiondegrees = [90,270,180,0,225,315,135,45];
subject = {'S01','S02','S03','S04','S05','S06'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
loction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
loc_num = [5,6,8,7,2,4,3,1];
edges= 0 : 22.5 : 360 ; 
time=[0,100];
% VU=[];
% VL=[];
% HL=[];
% HR=[];
% LL=[];
% LR=[];
% UL=[];
% UR=[];
total_d_loc = zeros(8,8,6);
for k =1 : 6 %subject
    VU=[];
    VL=[];
    HL=[];
    HR=[];
    LL=[];
    LR=[];
    UL=[];
    UR=[];
    count_edges = zeros(8,16,2);
    for i = 1 : 2 %condition 
        for j = 1 : 8 % direction
            main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{k}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
%             sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
            cd(MATpath)
    %         cd(fullfile(main_folder, 'eyedata','MATs'));
    %         ms_name = dir(sprintf('%s.mat', direction{j})).name;
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction{j}) );
            load(ms_path);
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
%             aaa=ls;
%             bb=aaa(strmatch(direction{j},aaa),:);
% %             b=bb(strmatch('tab_new_outside_blink',bb),:);
%             tab_path=fullfile(MATpath,bb(11,:));
%             load(tab_path);
            AM=sqrt(MS_TEMP(:,6).^2+MS_TEMP(:,7).^2); %add total amplitude
            MS_TEMP=[MS_TEMP,AM];
            MS_TEMP=calculate_degree(MS_TEMP);
            MS_TEMP(MS_TEMP(:,12)<0,12)=MS_TEMP(MS_TEMP(:,12)<0,12)+360;
            MS_TEMP(MS_TEMP(:,13)<0,13)=MS_TEMP(MS_TEMP(:,13)<0,13)+360;
            %find MS trial
            lower=(1300-time(1))*samplingRateData/1000;
            higher=(1800+time(2))*samplingRateData/1000;
            ms_start=MS_TEMP(:,8);
            ms_end=MS_TEMP(:,9);
            ms_trial=MS_TEMP(:,10);
            a= lower<ms_start & higher > ms_start;
            trial_w1=ms_trial(a);
            b= lower<ms_end & higher > ms_end;
            trial_w2=ms_trial(b);
            trial_w=unique([trial_w1;trial_w2]); 
            for ii = 1 :  length(trial_w)
                loc = tab(trial_w(ii),9);
                switch loc
                    case 1
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      UR = [UR;ind];
                    case 2
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      LL = [LL;ind];      
                    case 3
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      UL = [UL;ind];
                    case 4
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      LR = [LR;ind];  
                    case 5
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      VU = [VU;ind];
                    case 6
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      VL = [VL;ind];      
                    case 7
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      HR = [HR;ind];
                    case 8
                      ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                      HL = [HL;ind];  
                end
            end
%             cc=[];
%             for ii = 1 :  length(trial_w)
%                 ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
%                 cc=[cc;ind];
%             end
%             kk = histcounts(cc,edges);
%             count_edges(j,:,i)=kk;    
%             if i == 2
%                 if j < 5
%                     for ii = 1 :  length(trial_w)
%                         ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
%                         ind = ind - directiondegrees(j);
%                         con1 = [con1;ind];
%                     end
%                 else
%                     for ii = 1 :  length(trial_w)
%                         ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
%                         ind = ind - directiondegrees(j);
%                         con2 = [con2;ind];
%                     end
%                 end
%             else
%                 if j < 5
%                     for ii = 1 :  length(trial_w)
%                         ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
%                         ind = ind - directiondegrees(j);
%                         con3 = [con3;ind];
%                     end
%                 else
%                     for ii = 1 :  length(trial_w)
%                         ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
%                         ind = ind - directiondegrees(j);
%                         con4 = [con4;ind];
%                     end
%                 end 
%             end
%             %transfrom >180 to negative
%             con1(con1>180)=con1(con1>180)-360;
%             con2(con2>180)=con2(con2>180)-360;
%             con3(con3>180)=con3(con3>180)-360;
%             con4(con4>180)=con4(con4>180)-360;
%             
        end
    end
    total_d = zeros(8,8);
    VU_num=histcounts(VU,edges);
    VL_num=histcounts(VL,edges);
    HL_num=histcounts(HL,edges);
    HR_num=histcounts(HR,edges);
    LL_num=histcounts(LL,edges);
    LR_num=histcounts(LR,edges);
    UL_num=histcounts(UL,edges);
    UR_num=histcounts(UR,edges);
    aa1 = [VU_num;VL_num;HL_num;HR_num;LL_num;LR_num;UL_num;UR_num];
    for o = 1 : 7 
        total_d(:,o+1)=aa1(:,2*o)+aa1(:,2*o+1);
    end    
    total_d(:,1) = aa1(:,1)+aa1(:,16);    
    edges_plot=-0.125*pi:0.25*pi:1.875*pi;
    total_d_loc(:,:,k)=total_d;
 
    for pp = 1 : 8
        title_name = [subject{k},'_ ',direction{pp},' direction'];
        figpath='C:\Users\86186\Desktop\fig';
        name = [figpath,'/',subject{k},'_',direction{pp},'_dir_count_end_loc.png'];
        figure
        polarhistogram('BinEdges',edges_plot,'BinCounts',total_d(pp,:),'FaceAlpha',.1)
        
        title(title_name)
%        saveas(gca,name)
    end

    
%     figure
%     histogram(con1,'BinEdges',edges)
%     title('Cardinal directions ON meridians')
%     saveas(gca, sprintf('%s/degree1.png',MATpath))
% 
%     figure
%     histogram(con2,'BinEdges',edges)
%     title('Oblique directions OFF meridians')
%     saveas(gca, sprintf('%s/degree2.png',MATpath)) 
% 
%     figure
%     histogram(con3,'BinEdges',edges)
%     title('Cardinal directions OFF meridians')
%     saveas(gca, sprintf('%s/degree3.png',MATpath)) 
% 
%     figure
%     histogram(con4,'BinEdges',edges)
%     title('Oblique directions ON meridians')
%     saveas(gca, sprintf('%s/degree4.png',MATpath)) 
end
figpath='C:\Users\86186\Desktop\fig';
filename = [figpath,'/','location.mat'];
save(filename,'total_d_loc')

