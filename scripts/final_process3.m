clc; clear all;
homedir = '/Users/rania/Downloads/MS_Project';
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
sub_fir6 = {'S01','','S02','','S03','','S04','','S05','','S06'};
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
time=[0,100];
color = {[51, 34, 136],[17, 119, 51],[68, 170, 153],[136, 204, 238],[221, 204, 119],[204, 102, 119],[170, 68, 153],[136, 34 85]};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
%color = {'-r','-y','-g','-c','-b','-k'};
x = (1 : 2500) - 1300;
figure

allsubs_oblique_hard = []; vectorsubs_oblique_hard = []; % RE
allsubs_cardinal_hard = []; vectorsubs_cardinal_hard = []; 
allsubs_oblique_easy = [];  vectorsubs_oblique_easy = []; 
allsubs_cardinal_easy = []; vectorsubs_cardinal_easy = []; 

for kk = 1:8 % : 8
    
    if kk == 5 %|| kk == 5
        continue
    end
    
    sub_d = [];
    for ti= 1:2 % tilt easy / hard
        for i = 1 : 2
            if i == 1 & kk > 6
                continue
            end
            for j = 1 : 8 %8

                    if ti==1
                        tiltlevel = [0.5, 1, 2]; %, 1];
                    elseif ti==2
                        tiltlevel = [4, 8]; %[4, 8]; %, 8];
                    end

                    main_folder = fullfile(homedir,'Data_DI_wEYE', subject{kk}, ...
                         'RawData', condition{i}, 'Block1');
                    cd(fullfile(main_folder, 'eyedata'));
                    edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
                    edf_path = fullfile(main_folder,'eyedata',edf_name);
                    msg_filepath=replace(edf_path,'edf','msg');
                    samplingRateData=findSamplingRate(msg_filepath);
        %            sample=[sample,samplingRateData];
                    MATpath = fullfile(main_folder, 'eyedata','MATs');

                    ms_path= fullfile(MATpath,sprintf('%s.mat', direction_name{j}) );
                    load(ms_path);

                    tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
                    load(tab_path)

                    timepoint = countfigure_new(MS_TEMP,tab,samplingRateData, tiltlevel);
                    s_mean = nanmean(timepoint,1)*60; % convert to hz
                    sub_d = [sub_d;s_mean(1000:2194)];

                    if ti ==1
                        if j>4
                            allsubs_oblique_hard = [allsubs_oblique_hard; timepoint];
                        else
                            allsubs_cardinal_hard = [allsubs_cardinal_hard; timepoint];
                        end
                    elseif ti==2
                        if j>4
                            allsubs_oblique_easy = [allsubs_oblique_easy; timepoint];
                        else
                            allsubs_cardinal_easy = [allsubs_cardinal_easy; timepoint];
                        end
                    end


            end
        end
        hold on
        %a = shadedErrorBar(x(1000:2194), mean(sub_d,1),std(sub_d,0,1)/sqrt(size(sub_d,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
        %hold on
        if ti==1 % hard
            b = plot(x(1000:2194), mean(sub_d,1),'-','LineWidth',1);
            %hold on
            %obl_mean_hard_sem = nanstd(sub_d)./sqrt(size(sub_d,1));
            %a = shadedErrorBar(x(1000:2194), mean(sub_d,1),std(sub_d,0,1)/sqrt(size(sub_d,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
            %hold off
        else % easy
            b = plot(x(1000:2194), mean(sub_d,1),'-','LineWidth',3);
            %hold on
            %a = shadedErrorBar(x(1000:2194), mean(sub_d,1),std(sub_d,0,1)/sqrt(size(sub_d,1)), 'lineprops', {'-','Color',color{kk}/255},'transparent',1,'patchSaturation',0.2);
            %hold off
        end
         
        %a.Color = color{kk}/255;
        b.Color = color{kk}/255;
    
        % RE Added
        try
            temp_card = [sub_d(1:4,:);sub_d(9:12,:)]; 
            temp_obl = [sub_d(5:8,:);sub_d(13:16,:)]; 
        catch
            temp_card = [sub_d(1:4,:)]; 
            temp_obl = [sub_d(5:8,:)]; 
        end
        
       if ti ==1
            vectorsubs_oblique_hard = [vectorsubs_oblique_hard; mean(temp_obl,1)];
            vectorsubs_cardinal_hard = [vectorsubs_cardinal_hard; mean(temp_card,1)];
            %vectorsubs_oblique_hard = [vectorsubs_oblique_hard; temp_obl];
            %vectorsubs_cardinal_hard = [vectorsubs_cardinal_hard; temp_card];
        elseif ti==2
            vectorsubs_oblique_easy = [vectorsubs_oblique_easy; mean(temp_obl,1)];
            vectorsubs_cardinal_easy = [vectorsubs_cardinal_easy; mean(temp_card,1)];
            %vectorsubs_oblique_easy = [vectorsubs_oblique_easy; temp_obl];
            %vectorsubs_cardinal_easy = [vectorsubs_cardinal_easy; temp_card];
        end
        
    end % RE
        
end

% for hh = 1 : 6    
%     hold on
%     plot(x, sub_d(hh,:),color{hh},'LineWidth',1)
%     hold on
% end
xline(0,'-','Stimuli on')
hold on 
xline(500,'-','Stimuli off')
ylim([0,6])
xlim([-300,894])
%legend(sub_fir6)
ylabel('Microsaccade Rate (hz)')
xlabel('Time (ms)')
%title('Time series')

%%
figure
cond1 = [vectorsubs_cardinal_easy] ; %; vectorsubs_cardinal_hard];
cond2 = [vectorsubs_oblique_hard] ; 

windowSize = 10; %50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;

easy = filter(b,a,cond1); %imgaussfilt(cond1,3);
hard = filter(b,a,cond2); %imgaussfilt(cond2,3);
b_SEM = std(cond1-cond2)/sqrt(size(easy,1));

extremes = easy - hard;
%extremes = imgaussfilt(extremes,5);
yvec = mean(extremes,1); % just so +- SEM works quick
b = plot(x(1000:2194), yvec,'k-','LineWidth',2);
hold on
%b_SEM = imgaussfilt(b_SEM,1);
% plot(x(1000:2194), yvec+b_SEM,'r-','LineWidth',1);
% hold on
% plot(x(1000:2194), yvec-b_SEM,'r-','LineWidth',1);
hold on
a = shadedErrorBar(x(1000:2194), yvec, b_SEM, 'lineprops', {'-','Color',[0 0 0]/255},'transparent',1,'patchSaturation',0.2);
hold on
xline(0, 'LineWidth',2)
hold on
xline(500, 'LineWidth',2)
hold on
yline(0, 'r--', 'LineWidth',2)
ax = gca;
ax.XLim = [-100 800];
ax.YLim = [-.3 .75];
ax.XTick = [0, 200, 400, 600, 800];
ax.XTickLabel = {'0', '200', '400', '600', '800'};
ax = gca;
%ax.FontSize = 20;
ax.LineWidth = 2;
f1 = gcf;
f1.Position = [274 527 740 78];
box off

%figpath='C:\Users\86186\Desktop\fig\new\time_series.png';
%title_file = [figpath,'/',subject{kk},'- within stimuli.png'];
%saveas(gca,figpath)

%%

% sem across subjects

card_mean_hard = nanmean(allsubs_cardinal_hard,1)*60;
obl_mean_hard = nanmean(allsubs_oblique_hard,1)*60;
card_mean_easy = nanmean(allsubs_cardinal_easy,1)*60;
obl_mean_easy = nanmean(allsubs_oblique_easy,1)*60;

% card_mean_hard_sm = imgaussfilt(card_mean_hard,3);
% card_mean_easy_sm = imgaussfilt(card_mean_easy,3);
% obl_mean_hard_sm = imgaussfilt(obl_mean_hard,3);
% obl_mean_easy_sm = imgaussfilt(obl_mean_easy,3);

windowSize = 10; %50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
card_mean_hard_sm = filter(b,a,card_mean_hard');
card_mean_easy_sm = filter(b,a,card_mean_easy');
obl_mean_hard_sm = filter(b,a,obl_mean_hard');
obl_mean_easy_sm = filter(b,a,obl_mean_easy');

figure
b = plot(obl_mean_hard_sm,':','Color', [175, 141, 195]/255, 'LineWidth',2);
hold on
% a = shadedErrorBar(x(1000:2194), obl_mean_hard_sm, obl_mean_hard_sem, 'lineprops', {'-','Color',[0 0 0]/255},'transparent',1,'patchSaturation',0.2);
hold on
b = plot(obl_mean_easy_sm,'-','Color', [175, 141, 195]/255, 'LineWidth',2);
hold on
b = plot(card_mean_hard_sm,':','Color', [127, 191, 123]/255, 'LineWidth',2);
hold on
b = plot(card_mean_easy_sm,'-','Color', [127, 191, 123]/255, 'LineWidth',2);
xlim([1000 2194])
hold on
xline(1300, 'LineWidth',2)
hold on
xline(1800, 'LineWidth',2)
ylim([0 3])
ax = gca;
ax.XTick = [1100, 1300, 1500, 1700, 1900, 2100];
ax.XTickLabel = {'-200', '0', '200', '400', '600', '800'};
ax.YTick = [0 1 2 3];
f1 = gcf;
f1.Position = [78 303 740 182]; %[16 579 811 183];
ax.XLim = [1200 2100];
ax = gca;
%ax.FontSize = 20;
ax.LineWidth = 2;
box off
