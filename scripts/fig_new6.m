clear all; clc;
%for stimulus direction plot
% directiondegrees = [45,225,135,315,90,270,0,180];
directiondegrees = [90,270,180,0,225,315,135,45];
subject = {'S01','S02','S03','S04','S05','S06'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
edges= 0 : 22.5 : 360 ; 
time=[0,100];
total_d_dir1 = zeros(8,8,6);
for k =1 : 1 %subject
    count_edges = zeros(8,16,2);
    count_acc = zeros(8,2,16,2); %direction * condtion * accur (16*2) 
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
            cc=[];
            trail_T = [];
            pre = [];
            for ii = 1 :  length(trial_w)
                ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                trail_num = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),10);
                cc=[cc;ind];
                trail_T = [trail_T;trail_num];
            end
            [kk,~,indx] = histcounts(cc,edges);
            accur = zeros(16,2);
            for osc = 1 : 16 %for 16 interval
                oli = unique(trail_T(indx==osc));
                Acc = [];
                for jia = 1 : length(oli)
                    Acc = [Acc; tab(tab(:,1)==oli(jia),15)];
                    
                end
                accur(osc,:) =  [sum(Acc);length(Acc)];
                pre = [pre; mean(Acc)];
            end
            count_acc(j,i,:,:)=accur  ;                     
            count_edges(j,:,i)=kk;                
        end
    end
    count_acc1 = squeeze(count_acc(:,1,:,:)+count_acc(:,2,:,:));
    count_acc2 = zeros(8,8,2);
    for o = 1 : 7 
        count_acc2(:,o+1,:)=count_acc1(:,2*o,:)+count_acc1(:,2*o+1,:);
    end
    count_acc2(:,1,:)=count_acc1(:,1,:)+count_acc1(:,16,:);
        
    count_acc3 = squeeze(count_acc2(:,:,1)./count_acc2(:,:,2));
%    total_d_acc = zeros(8,8);
    
    aa1= count_edges(:,:,1)+count_edges(:,:,2);
    total_d_dir = zeros(8,8);
    for o = 1 : 7 
        total_d_dir(:,o+1)=aa1(:,2*o)+aa1(:,2*o+1);
 %       total_d_acc(:,o+1)=count_acc2(:,2*o)+count_acc2(:,2*o+1);
    end
    total_d_dir(:,1) = aa1(:,1)+aa1(:,16);
    total_d_dir = total_d_dir./sum(total_d_dir,2);
    count_acc3 = count_acc3./sum(count_acc3);
    edges_plot=-0.125*pi:0.25*pi:1.875*pi;
    total_d_dir1(:,:,k)=total_d_dir;
    for pp = 1 : 8 
        title_name = [subject{k},'_ ',direction{pp},'_ dir/acc'];
        figpath='C:\Users\86186\Desktop\fig';
        name = [figpath,'/',subject{k},'_',direction{pp},'_dir_count_acc.png'];
        figure
        polarhistogram('BinEdges',edges_plot,'BinCounts',total_d_dir(pp,:),'FaceAlpha',.0,'LineWidth',1.5,'EdgeColor','r')
        hold on 
        polarhistogram('BinEdges',edges_plot,'BinCounts',count_acc3(pp,:),'FaceAlpha',.0,'LineWidth',1.5,'EdgeColor','b')
        legend('% of MS','Accuracy')
        title(title_name)
        saveas(gca,name)
    end
end
%filename = [figpath,'/','direction.mat'];
%save(filename,'total_d_dir1')

