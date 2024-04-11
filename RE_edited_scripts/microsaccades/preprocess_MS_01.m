clc; clear all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')

warning('off', 'MATLAB:interp1:NaNstrip');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
% Convert EDF and Preprocess eyedata 
% Read JSON file from parent directory for paths and params
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
datadir = jsonParams.datadir.Path; 
scriptsdir = jsonParams.scriptsdir.Path;
edf2asc = jsonParams.edf2ascdir.Path;

addpath(genpath(scriptsdir))

subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};

% Best to use in debug mode
plotON = 1; % this will create plots for each trial iteratively on the same figure handle.

% set a predetermined temporal length for event matrix
nSampleCutOff = 4000; %2500;  % this is in ms & samples (data will always be in 1000hz)

VFAC=[6,6]; % velocity threshold for X, Y to constitute as a microsaccade
MINDUR=6; % minimum duration (ms) to constitute as a microsaccade 
mergeInterval=15; % minimum duration (ms) to constitute as separate microsaccades; otherwise merged as one 
threshold_AM=[.05 1]; %[0.05,1]; % minimum and maximum amplitude (deg) to constitute as a microsaccade 

%%

tic
for ii = 2 : 2
    for jj = 2 : 2 %2 %2 % NOT 1:2 (for S07 / S08)
        for kk = 1 : 8 %8 %8 %8 %1 : 8 %8

            % 2 out of 112 eyelink files are corrupted (missing data).
            % Leave them out for the analysis
            if ii==7 && jj==2 && kk==2 % this eyelink file is corrupted (missing data)
                continue
            elseif ii==8 && jj==2 && kk==5 % this eyelink file is corrupted (missing data)
                continue
            end

            rawdata_allblocks = fullfile(datadir, subject{ii}, ...
                'RawData', condition{jj});
            processeddata_folder = fullfile(datadir, subject{ii}, ...
                'ProcessedData', condition{jj});
            processedeyeddata_folder = fullfile(datadir, subject{ii}, ...
                'ProcessedData', condition{jj}, 'eyedata','MATs');
            figures_folder = fullfile(datadir, subject{ii}, ...
                'ProcessedData', condition{jj}, 'eyedata','figures');
            microsaccdata_folder = fullfile(datadir, subject{ii}, ...
                'ProcessedData', condition{jj}, 'eyedata','microsaccades');


            % check for possibility of multiple blocks
            subdirectories = checkBlockn(rawdata_allblocks);
            numBlocks = length(subdirectories);

            for bi=1:numBlocks

                % for cases with repeated sessions, save both.
                if bi==1
                    blocklabel = '';
                else
                    blocklabel = num2str(bi);
                end

                rawdata_folder = fullfile(rawdata_allblocks, subdirectories{bi});
    
                csv_filepath = fullfile(rawdata_folder,sprintf('expRes%s_1dMotionAsym_Psychophysics_%s.csv', ...
                    subject{ii}, direction{kk}));
                scr_filepath = fullfile(rawdata_folder,sprintf('scr_file%s_1dMotionAsym_Psychophysics_%s.mat', ...
                    subject{ii}, direction{kk}));

                try
                    cd(fullfile(rawdata_folder, 'eyedata')); edf_name = dir(sprintf('*%s*.edf', direction{kk})).name;
                    edf_path = fullfile(rawdata_folder,'eyedata',edf_name);
                    edfExists = 1;
                catch
                    edfExists = 0;
                end

                if edfExists
                
                    [path_folder,~,~] = fileparts(edf_path);
        
                    % make processed data directories if does not exist
                    mkdir(processedeyeddata_folder);
        
                    % load the screen properties
                    load(scr_filepath)
                    distanceFromScreen = scr.dist;      % dist of monitor from eye (in cm)
                    screenWidthCm = scr.disp_sizeX/10;  % width of monitor (in cm from mm)
                    screenWidthPx = scr.scr_sizeX;      % number of pixels (horizontally)
                    screenHeightPx = scr.scr_sizeY;     % number of pixels (vertically)
                    screenCenter = [screenWidthPx/2 screenHeightPx/2];  % screen center (intial fixation position)
                    dvaPerPx = rad2deg(atan2(.5*screenWidthCm, distanceFromScreen)) / (.5*screenWidthPx);
        
                    % convert EDF (if not already converted
                    if ~isfile(replace(edf_path,'edf','msg'))
                        disp('Converting EDF.. please wait.')
                        % creates initial _Dat_stim, _Dat_all, _tab and _blink
                        [path_tab,path_blink,path_stim,path_dat, samplingRateData]=convert(edf_path, ...
                            edf2asc,processedeyeddata_folder);
                        % updates _tab with behavioral data per trial
                        converge(path_tab,csv_filepath);
                        % checks Dat_stim for eye position outside 1.5 deg
                        check_outside(path_stim, screenCenter, dvaPerPx);
                        % updates _tab file with boolean if eye went outside
                        % 1.5 deg during stim period
                        add_outside(path_tab,path_stim,samplingRateData);
                        % updates _tab file with boolean if eye blinked during
                        % stim period
                        add_blink(path_tab,path_blink, samplingRateData);
                        % updates _dat_all file with boolean if eye blinked 
                        omitblinks(path_dat,path_blink,0,samplingRateData);
                        disp('Complete.')
                    else
                        % if already converted, just get sampling rate
                        msg_filepath=replace(edf_path,'edf','msg');
                        samplingRateData=findSamplingRate(msg_filepath);
                        path_dat = fullfile(processedeyeddata_folder, replace(edf_name, '.edf', '_Dat_all.mat'));
                        path_tab = fullfile(processedeyeddata_folder, replace(edf_name, '.edf', '_tab.mat'));
                    end
    
                    load(path_tab)
                    load(path_dat)
        
                    %% Analyze microsaccades
        
                    % check tab/dat compatibility:
                    % in rare cases, eyelink file gets too large and stops
                    % recording data even though the message file exists
                    % any tab entry should be less than the last time stamp
                    % in the dat_all file
                    cc = Dat_all(end,1) < tab(:,8); 
                    outofRange_idx = find(cc == 1, 1, 'first');
                    if ~isempty(outofRange_idx)
                        if tab(end,8)-Dat_all(end,1) < 2 % sometimes it is off by 1 ms (just repeat the last point to include the trial
                            Dat_all(end+1,:) = Dat_all(end,1);
                        else
                            tab = tab(1:outofRange_idx-1, :); % do not save this file, better to have original
                            % this is confirmed for S03-2-4, S03-2-7,
                            % S07-2-2, S08-2-5
                            % (only ~ half trials are saved)
                            warning('Dat_all appears to be shorter than trial stamps. Check the data before assuming data loss.')
                        end
                    end

                    % # of trials
                    [nTrials, ~] = size(tab);
                    numBlink=[]; MS=[];
        
                    eventMatrix = nan(nTrials, nSampleCutOff);
                    eyetraceMatrix = nan(nTrials, nSampleCutOff,2);
        
                    % these are not used
                    col_dva=nan(size(Dat_all,1),2);
                    velo=nan(size(Dat_all,1),2);
        
                    msgFileSamplingRateData = samplingRateData; % log this value to initiate every iteration of the for loop
        
                    % loop through # of trials
                    for i = 1 : nTrials
        
%                         % this was here before the if statement was added
%                         % inside the top of the for loop:
%                         if tab(i,8) == 0 % very rare cases, the eyetracker timestamp gliches to 0
%                             sprintf('Dropping trial %i faulty trial', i)
%                             sprintf('TRIAL_END appears right after the TRIAL_START: %i', tab(i,2))
%                             continue
%                         end
        
                        % get Dat for trial i (4-5 cols)
                        order = tab(i,2) < Dat_all(:,1) & tab(i,8) > Dat_all(:,1);
                        trial= Dat_all(order,:);
        
                        samplingRateData = msgFileSamplingRateData;
        
                        % upsample or downsample in cases of different
                        % sampling rate (e.g., 1000hz vs 2000hz):
                        % this will also change samplingRateData
                        [trial,samplingRateData] = check4interp(trial,samplingRateData);
                       
                        % make labels for each eye trace segment in between blinks
                        new_trial=segmentnonBlinks2(trial);
        
%                         % check for rare eyetracker error (I think I already do
%                         % this above)
%                         if ~isempty(new_trial)
%                             if (range(new_trial(:,2)) == 0) && (range(new_trial(:,3)) == 0) && (range(new_trial(:,4)) == 0)
%                                 new_trial(:, 5) = NaN;
%                                 new_trial(:, 6) = NaN;
%                             end
%                         end
        
                        % count the number of segments between blinks
                        k=new_trial(:,6);
                        num_seg=unique(k(~isnan(new_trial(:,5))));
                        numBlink=[numBlink,size(num_seg,1)];
        
                        % fill event matrix (cut off if longer than )
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
                        if isempty(MS)
                            continue
                        end
        
                        % add label 3 for saccades in event matrix
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
        
                        % add label 2 for micrpsaccades in event matrix
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
        
                        if plotON
                            % overlay all MS for this trial
                            ms2plot = MS_fil(MS_fil(:,10)==i,:);
                            ms2plot = ms2plot(ms2plot(:,2)<=nSampleCutOff,:); % remove MS outside of 2500 range
                            checkQuality(screenCenter, dvaPerPx, new_trial, ms2plot, i, tab)
                        end
                        
                    end 
        
                    %%
        
                    % compute rate and event data 
                    EVENTS = eventMatrix(:,1:nSampleCutOff);
                    EVENTS(EVENTS==0) = NaN; %blinks
                    EVENTS(EVENTS==3) = NaN; %saccades
        
                    EVENTS = EVENTS-1; % 2 to 1
        
                    rate = nansum(EVENTS,1)./sum(~isnan(EVENTS),1);
        
                    mkdir(figures_folder)
                    mkdir(microsaccdata_folder)
        
                    figure
                    imagesc(EVENTS)
                    saveas(gca, sprintf('%s/MS_events_%s%s.png',figures_folder, direction{kk},blocklabel))
                    figure
                    plot(rate)
                    hold on
                    xline(1300, 'r')
                    hold on 
                    xline(1800, 'r')
                    saveas(gca, sprintf('%s/MS_rate_%s%s.png',figures_folder, direction{kk},blocklabel))
        
                    eventFilename = strcat(subject{ii},'_Switched_',direction{kk}, blocklabel, '_events.mat');
                    save(fullfile(microsaccdata_folder,eventFilename), 'EVENTS');

                    eventFilename = strcat(subject{ii},'_Switched_',direction{kk}, blocklabel, '_allevents.mat');
                    save(fullfile(microsaccdata_folder,eventFilename), 'eventMatrix');
        
                    rateFilename = strcat(subject{ii},'_Switched_',direction{kk}, blocklabel, '_rate.mat');
                    save(fullfile(microsaccdata_folder,rateFilename), 'rate');
                    
                    
                    %% MS Characteristics
        
                    % leave out microsaccades that overlap with the end of
                    % nSampleCutOff (e.g., 2500 ms post trial onset)
                    MS_TEMP=getms(MS_fil,samplingRateData,nSampleCutOff);
        
                    AM=sqrt(MS_TEMP(:,6).^2+MS_TEMP(:,7).^2); %add total amplitude
                    MS_TEMP=[MS_TEMP,AM];
        
                    % calculate direction (0-360 polar angle)
                    MS_TEMP=calculate_direction_RE(MS_TEMP);
        
                    edge=0.125:0.25:2.125;
                    edge=edge-1;
                    edge1=edge*pi;
                    edge2=edge*180;
        
                    %% MS direction
        
                    % selection of timepoints to analyze microsaccade direction
                    stimpost_ms = MS_TEMP(MS_TEMP(:,9)>1.3*samplingRateData,:); % end point after stimulus pres
                    stimpost_ms = stimpost_ms(stimpost_ms(:,8)<1.8*samplingRateData,:); % start point before  stim offset
        
                    rho_degs = linspace(0,359.99,15);
        
                    figure
                    polarhistogram(deg2rad(stimpost_ms(:,12)),deg2rad(rho_degs))
                    saveas(gca, sprintf('%s/MS_dir_amp_%s%s.png',figures_folder, direction{kk},blocklabel))
                    figure
                    polarhistogram(deg2rad(stimpost_ms(:,13)),deg2rad(rho_degs))
                    saveas(gca, sprintf('%s/MS_dir_component_%s%s.png',figures_folder, direction{kk},blocklabel))
        
                    % save matrix of MS data
                    save(sprintf('%s/%s%s_microsaccadeMatrix.mat',microsaccdata_folder, direction{kk},blocklabel),'MS_TEMP');
        
                    dirRespCongruent = 0;
                    dirRespIncongruent = 0;
                    dirRespCongruent_incTrials = 0;
                    dirRespIncongruent_incTrials = 0;
                    dirRespCongruent_cTrials = 0;
                    dirRespIncongruent_cTrials = 0;

                    % just added
                    % log how the stimpost_ms direction corresponds with behavior:
                    [numStimMS, ~] = size(stimpost_ms);
                    for ss=1:numStimMS
                        trialID = stimpost_ms(ss,10);

                        % standardDir
                        directiondegrees = {45, 225, 135, 315, 90, 270, 0, 180};
                        standardDir = directiondegrees{tab(trialID,10)};

                        % correct for negatives (e.g., -358 is -2?
                        if tab(trialID,12)<0
                            curDir = (360-(tab(trialID,12)))*-1;
                        else
                            curDir = tab(trialID,12);
                        end

                        % TARGET WITH RESPECT TO THE STANDARD

                        % Convert angles to radians
                        stim_rad = deg2rad(standardDir);
                        target_rad = deg2rad(curDir);

                        % Calculate angular distance from stim_rad to ms_rad
                        angular_targetdistance = target_rad - stim_rad;

                        % Normalize the angular distance to be between -pi and pi
                        angular_targetdistance = mod(angular_targetdistance + pi, 2*pi) - pi;
                        if angular_targetdistance > 0
                            target_clockwise = 0;
                            if tab(trialID, 14) == 1 % if correct
                                respC = 0;
                            else 
                                respC = 1;
                            end
                        elseif angular_targetdistance < 0
                            target_clockwise = 1;
                            if tab(trialID, 14) == 1 % if correct
                                respC = 1;
                            else 
                                respC = 0;
                            end
                        end

                        % MS WITH RESPECT TO THE STANDARD (THROUGH ORIGIN)
                        ms_rad = deg2rad(stimpost_ms(ss, 12)); % 12 or 13?
                                                
%                         % check is diff is within 90 deg
%                         angular_distance = rad2deg(abs(ms_rad - stim_rad)); 
%                         
%                         if angular_distance <= 90
%                             if tab(trialID, 14) == 1  % MS dir is within 180
%                                 dirRespCongruent_cTrials = dirRespCongruent_cTrials+1;
%                             elseif tab(trialID, 14) == 0
%                                 dirRespCongruent_incTrials = dirRespCongruent_incTrials+1;
%                             end
%                         else 
%                             if tab(trialID, 14) == 1  % MS dir is within 180
%                                 dirRespIncongruent_cTrials = dirRespIncongruent_cTrials+1;
%                             elseif tab(trialID, 14) == 0
%                                 dirRespIncongruent_incTrials = dirRespIncongruent_incTrials+1;
%                             end
%                         end

                        angular_distance = ms_rad - stim_rad; % check whether point is c/cc
                        angular_distance = mod(angular_distance + pi, 2*pi) - pi;

                        % Determine if the MS direction is clockwise or counterclockwise
                        if angular_distance > 0 && ~respC ... % MS dir and response are counterclockwise
                               || angular_distance < 0 && respC % MS dir and response are clockwise 
                            dirRespCongruent = dirRespCongruent+1;
                            if tab(trialID, 14) == 0 % separately count the wrong trials
                                dirRespCongruent_incTrials = dirRespCongruent_incTrials+1;
                                %disp('WRONG AND CONGRUENT!')
                            elseif tab(trialID, 14) == 1
                                dirRespCongruent_cTrials = dirRespCongruent_cTrials+1;
                                %disp('CORRECT AND CONGRUENT!')
                            end
                        else % they are incongruent
                            dirRespIncongruent = dirRespIncongruent+1;
                            if tab(trialID, 14) == 0 % separately count the wrong trials
                                dirRespIncongruent_incTrials = dirRespIncongruent_incTrials+1;
                                %disp('WRONG AND INCONGRUENT!')
                            elseif tab(trialID, 14) == 1
                                dirRespIncongruent_cTrials = dirRespIncongruent_cTrials+1;
                                %disp('CORRECT AND INCONGRUENT!')
                            end
                        end
                        
                    end

                    %disp(sprintf('Percent of MS Congruent with Response: %d', dirRespCongruent/(dirRespCongruent+dirRespIncongruent)))
                    disp(sprintf('Percent of MS-Incorrect Trials with Incongruent MS: %.2f', dirRespIncongruent_incTrials/(dirRespIncongruent_incTrials+dirRespCongruent_incTrials)))
                    disp(sprintf('Percent of MS-Correct Trials with Incongruent MS: %.2f', dirRespIncongruent_cTrials/(dirRespIncongruent_cTrials+dirRespCongruent_cTrials)))
                    disp(sprintf('Percent of MS-Congruent Trials that are Correct: %.2f', dirRespCongruent_cTrials/(dirRespIncongruent_incTrials+dirRespCongruent_cTrials)))
                    disp(sprintf('Percent of MS-Incongruent Trials that are Correct: %.2f', dirRespIncongruent_cTrials/(dirRespIncongruent_incTrials+dirRespIncongruent_cTrials)))
                    '~~~~~~~~~~~~~~~~~'

                    %% Summary stats of microsaccades for the session
                    summary_trial=count_ms(MS_TEMP,tab,microsaccdata_folder, direction{kk});
                    save(sprintf('%s/%s%s_summary.mat',microsaccdata_folder,direction{kk},blocklabel),'summary_trial');
    
                    close all;
                end
            end
        end
    end
end
toc