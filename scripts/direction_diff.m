clear all; clc;
%for direction plot
% directiondegrees = [45,225,135,315,90,270,0,180];
directiondegrees = [90,270,180,0,225,315,135,45];
subject = {'S01','S02','S03','S04','S05','S06'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
edges= -180 : 45 : 180 ; 
time=[0,100];
con1 = [];
con2 = [];
con3 = [];
con4 = [];
for k =1 : 6 %subject
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
            if i == 2
                if j < 5
                    for ii = 1 :  length(trial_w)
                        ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                        ind = ind - directiondegrees(j);
                        con1 = [con1;ind];
                    end
                else
                    for ii = 1 :  length(trial_w)
                        ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                        ind = ind - directiondegrees(j);
                        con2 = [con2;ind];
                    end
                end
            else
                if j < 5
                    for ii = 1 :  length(trial_w)
                        ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                        ind = ind - directiondegrees(j);
                        con3 = [con3;ind];
                    end
                else
                    for ii = 1 :  length(trial_w)
                        ind = MS_TEMP(trial_w(ii) == MS_TEMP(:,10),12);
                        ind = ind - directiondegrees(j);
                        con4 = [con4;ind];
                    end
                end 
            end
            %transfrom >180 to negative
            con1(con1>180)=con1(con1>180)-360;
            con2(con2>180)=con2(con2>180)-360;
            con3(con3>180)=con3(con3>180)-360;
            con4(con4>180)=con4(con4>180)-360;
            
        end
    end
    figure
    histogram(con1,'BinEdges',edges)
    title('Cardinal directions ON meridians')
    saveas(gca, sprintf('%s/degree1.png',MATpath))

    figure
    histogram(con2,'BinEdges',edges)
    title('Oblique directions OFF meridians')
    saveas(gca, sprintf('%s/degree2.png',MATpath)) 

    figure
    histogram(con3,'BinEdges',edges)
    title('Cardinal directions OFF meridians')
    saveas(gca, sprintf('%s/degree3.png',MATpath)) 

    figure
    histogram(con4,'BinEdges',edges)
    title('Oblique directions ON meridians')
    saveas(gca, sprintf('%s/degree4.png',MATpath)) 
end


