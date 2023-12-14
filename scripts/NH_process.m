clc; clear all;

subject = 'S03'; condition = 'Full_distance_non_radialtangential'; direction = 'VU';
%home_folder = '/Users/rania/Desktop/RadialBias_pilot1/Microsaccade_Analysis';
edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';

main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject, ...
    'RawData', condition, 'Block1');
csv_filepath = fullfile(main_folder,sprintf('expRes%s_1dMotionAsym_Psychophysics_%s.csv', ...
    subject, direction));
scr_filepath = fullfile(main_folder,sprintf('scr_file%s_1dMotionAsym_Psychophysics_%s.mat', ...
    subject, direction));
cd(fullfile(main_folder, 'eyedata')); edf_name = dir(sprintf('*%s*.edf', direction)).name;
edf_path = fullfile(main_folder,'eyedata',edf_name);
[path_folder,~,~] = fileparts(edf_path);

if ~isfile(replace(edf_path,'edf','msg'))
    [path_tab,path_stim,path_blink]=convert(edf_path,edf2asc);%create raw tab file and stim file and blink file
    path_newtab=converge(path_tab,csv_filepath);
    outside_path=check_outside(path_stim);
    outside_tab=add_outside(path_newtab,outside_path);
    final_tab=add_blink(outside_tab,path_blink);
end

msg_filepath=replace(edf_path,'edf','msg');
samplingRateData=findSamplingRate(msg_filepath);
MATpath = fullfile(main_folder, 'eyedata','MATs');

%path_blink='F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Experimental_SetUp\Data\NH\Full_distance_non_radialtangential\Block1\eyedata\MATs\VU120209_blink.mat';
path_blink = fullfile(MATpath, replace(edf_name, '.edf', '_blink.mat'));
%data_path='F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Experimental_SetUp\Data\NH\Full_distance_non_radialtangential\Block1\eyedata\MATs\VU120209_Dat_all.mat';
data_path = fullfile(MATpath, replace(edf_name, '.edf', '_Dat_all.mat'));
new_data=omitblinks(data_path,path_blink,0,samplingRateData);
%tab_path='F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Experimental_SetUp\Data\NH\Full_distance_non_radialtangential\Block1\eyedata\MATs\VU120209_tab_new_outside_blink.mat';
tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));

load(tab_path)
load(new_data)
load(scr_filepath)

distanceFromScreen = scr.dist;      % dist of monitor from eye (in cm)
screenWidthCm = scr.disp_sizeX/10;  % width of monitor (in cm from mm)
screenWidthPx = scr.scr_sizeX;      % number of pixels (horizontally)
screenHeightPx = scr.scr_sizeY;     % number of pixels (vertically)
screenCenter = [screenWidthPx/2 screenHeightPx/2];  % screen center (intial fixation position)
%dvaPerPx = atan2(1,distanceFromScreen)*180/pi/screenWidthPx * screenWidthCm; % degrees per pixel
% cm of my display; X=30.5 x Y=22.9
%h = 30.5; d=60; r = 1152;
dvaPerPx = rad2deg(atan2(.5*screenWidthCm, distanceFromScreen)) / (.5*screenWidthPx);

%%
nTrials = 800;
nSampleCutOff = 2.5*samplingRateData;
eventMatrix = nan(nTrials, nSampleCutOff);
eyetraceMatrix = nan(nTrials, nSampleCutOff,2);

col_dva=nan(size(Dat_all,1),2);
velo=nan(size(Dat_all,1),2);

VFAC=[6,6]; 
MINDUR=6;
mergeInterval=15;
threshold_AM=[0.05,1];

numBlink=[];
MS=[];
% figure
for i = 1 : 800
    order = tab(i,2) < Dat_all(:,1) & tab(i,8) > Dat_all(:,1);
    trial= Dat_all(order,:);
    new_trial=segmentnonBlinks2(trial);
    
    % check for rare eyetracker error
    if ~isempty(new_trial)
        if (range(new_trial(:,2)) == 0) && (range(new_trial(:,3)) == 0) && (range(new_trial(:,4)) == 0)
            new_trial(:, 5) = NaN;
            new_trial(:, 6) = NaN;
        end
    end

    k=new_trial(:,6);
    num_seg=unique(k(~isnan(new_trial(:,5))));
    numBlink=[numBlink,size(num_seg,1)];
    
    % fill event matrix
    temp = new_trial(:,5)';
    if length(temp) > nSampleCutOff
        temp = temp(1:nSampleCutOff);
    end
    eventMatrix(i, 1:length(temp)) = ~isnan(temp); % blinks (nan to 0)

    for j = 1 : size(num_seg,1)
        ord=k==num_seg(j);
        o=find(ord);
        trial_seg=new_trial(ord,:);
        if sum(isnan(trial_seg(:,5)))==0
            x=dvaPerPx*(trial_seg(:,2)-screenCenter(1));
            y=dvaPerPx*(trial_seg(:,3)-screenCenter(2));
            if isnan(y(end))
                y(end) = y(end-1);
            end
            if size(x,1) > 105
                x_fil=filtfilt(fir1(35,0.05),1,x);
                y_fil=filtfilt(fir1(35,0.05),1,y);
%             col_dva(o,1)=x_fil;
%             col_dva(o,2)=y_fil;
                d=[x_fil,y_fil];
                v = computevelocity(d,samplingRateData); 
            
            %velo(o,:)=v;
                %[msac, radius] = microsaccMerge_absolute(d,v,3, VFAC,MINDUR,mergeInterval,samplingRateData);
                [msac, radius] = microsaccMerge(d,v,VFAC,MINDUR,mergeInterval,samplingRateData);
                if isempty(msac)==0
                    num_seg_start=msac(:,1);
                    num_seg_end=msac(:,2);
                    trail_start=num_seg_start+o(1)-1;
                    trail_end=num_seg_end+o(1)-1;
                    num_trial = ones(size(num_seg_start))*i;
                    msac=[msac,trail_start,trail_end,num_trial];
                    MS=[MS;msac];
                    MS_fil=filterAM(MS,threshold_AM);
                end
            end
        end
%         col_dva(o,1)=x_fil;
%         col_dva(o,2)=y_fil;
%         Dat_all=[Dat_all,col_dva,velo];
    end
    
    % add S for specific trial (as saccades)
    trial_saccades = MS(MS(:,10)==i, :); [segs, ~] = size(trial_saccades); 
    if segs~=0
        for ti=1:segs
            saccadeStart = trial_saccades(ti, 8);
            saccadeEnd = trial_saccades(ti, 9);
            if saccadeEnd > nSampleCutOff
                saccadeEnd = nSampleCutOff;
            end
            eventMatrix(i, saccadeStart:saccadeEnd) = 3;
        end
    end
    
    % add MS for specific trial
    trial_saccades = MS_fil(MS_fil(:,10)==i, :); [segs, ~] = size(trial_saccades);
    if segs~=0
        for ti=1:segs
            saccadeStart = trial_saccades(ti, 8);
            saccadeEnd = trial_saccades(ti, 9);
            if saccadeEnd > nSampleCutOff
                saccadeEnd = nSampleCutOff;
            end
            eventMatrix(i, saccadeStart:saccadeEnd) = 2;
        end
    end
    
%     if size(x,1) > 105
%         x_fulltrial=dvaPerPx*(new_trial(:,2)-screenCenter(1));
%         y_fulltrial=dvaPerPx*(new_trial(:,3)-screenCenter(2));
%         %plot(filtfilt(fir1(35,0.05),1,x_fulltrial), 'b')
%         hold on
%         plot(filtfilt(fir1(35,0.05),1,y_fulltrial), 'r')
%         ylim([-2 2])
%     end
%     hold on
end 

%%
EVENTS = eventMatrix(:,1:nSampleCutOff);
EVENTS(EVENTS==0) = NaN; %blinks
EVENTS(EVENTS==3) = NaN; %saccades

EVENTS = EVENTS-1; % 0 to 1
event_name = [subject{ii},'_',c{jj},'_',direction{kk},'_events.mat'];
save(event_name,'EVENTS')
rate_name = [subject{ii},'_',c{jj},'_',direction{kk},'_rate.mat'];
rate = nansum(EVENTS,1)./sum(~isnan(EVENTS),1);
save(rate_name,'rate')

% figure
% imagesc(EVENTS)
% saveas(gca, sprintf('%s/MS_events_%s.png',MATpath, direction))
% figure
% plot(rate)
% saveas(gca, sprintf('%s/MS_rate_%s.png',MATpath, direction))
% %%
% MS_TEMP=getms(MS_fil,samplingRateData);
% save(sprintf('%s/%s.mat',MATpath, direction),'MS_TEMP');
% %%
% 
% time = samplingRateData/1000; % added these 3 lines as work around
% fil=MS_fil(:,9)/time < nSampleCutOff;
% ms_new=MS(fil,:);
% 
% %ms_new=drawfigure(MS_fil,samplingRateData,tab,data_path);
% %countfigure(ms_new,tab,samplingRateData,data_path)
% 
% AM=sqrt(ms_new(:,6).^2+ms_new(:,7).^2); %add total amplitude
% ms_new=[ms_new,AM];
% ms_new=calculate_degree(ms_new);
% 
% edge=0.125:0.25:2.125;
% edge=edge-1;
% edge1=edge*pi;
% edge2=edge*180;
% 
% %%
% 
% stimpost_ms = ms_new(ms_new(:,9)>1.3*samplingRateData,:); % end point after stimulus pres
% stimpost_ms = stimpost_ms(stimpost_ms(:,8)<1.8*samplingRateData,:); % start point before  stim offset
% 
% rho_degs = [-22.5 22.5 67.5 112.5 157.5 202.5 247.5 292.5 337.49999];
% figure
% polarhistogram(deg2rad(stimpost_ms(:,12)),deg2rad(rho_degs))
% saveas(gca, sprintf('%s/MS_dir_amp_%s.png',MATpath, direction))
% figure
% polarhistogram(deg2rad(stimpost_ms(:,13)),deg2rad(rho_degs))
% saveas(gca, sprintf('%s/MS_dir_component_%s.png',MATpath, direction))
% 
% MS_dirs=stimpost_ms(:,12:13);
% save(sprintf('%s/MS_dirs_%s.mat',MATpath, direction),'MS_dirs');
% %draw_his(ms_new,data_path)
% 
% %%
% summary_trial=count_ms(ms_new,tab,data_path);
% figure
% polarhistogram(ms_new(:,13)/180*pi,'BinEdges',edge1)
% 
% AM=sqrt(MS_HL(:,6).^2+MS_HL(:,7).^2); %add total amplitude
% MS_HL=[MS_HL,AM];
% de=atan(ms_new(:,6)./ms_new(:,7));%use amplitude 
% de2=atand(ms_new(:,6)./ms_new(:,7));
% ms_new=[ms_new,de,de2];
% 
% de=atan(ms_new(:,4)./ms_new(:,5));%use amplitude 
% de2=atand(ms_new(:,4)./ms_new(:,5));
% ms_new=[ms_new,de,de2];

% summary_trial=count_ms(ms_new,tab,edf_path);