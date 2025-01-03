% for pupil, retrieve values for current trial and next trial from ASC.
clc; clear all; close all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
%datadir = '/Users/rje257/Desktop/MS_Project/Data_DI_wEYE'; % override

%%
% PUPIL RESPONSE SHOULD BE DEFINED BY ASC PER TRIAL (trial_start to responseTime + 4500 ms)?

% EVENTS MATRIX (TIME ACROSS TRIAL): 
% (yes) Stimulus on (1300-1800ms), 
% (yes but not in tab file) saccade onset/latency?
% (yes) next trial stimulus on (fixation color change to black)
% (yes; in allevents.mat) [ENABLE: filter trials with ms; saccades; no ms during this period?]
% (yes; in allevents.mat) [ENABLE: filter out trials with SACCADE or EYE CLOSED?
% (yes; in allevents.mat - 0) [ENABLE: filter trials with BLINKS] - just
% interpolate linearlly

% 1 value per trial (tab file)
% Condition: tilt (-8 —> 8)
% Location (PA)
% Motion Direction (1-8)
% Session (1-16)
% response time, 
% response C/IC (feedback), 

direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};

%%

% for blink correction
pad_size = 150;
samplingRate = 1000; % hertz
stimStart = 1300;
stimEnd = 1800;
trialwise=1;
plotOn=0;

% this 0-centers the data based on this window
minWindowStart = 1000; minWindowEnd = 1800;
% this is used to calculate PSC
baselinePSCstart = 1000; baselinePSCend = 1300;

%%
for ii = 1 : 1 % CHANGE BACK TO 1:2 if running s04!!
    subj = subjects{ii};
    output = []; % every row is trial eventually all output is saved for this subject
    alltab = {};
    allfilteredPupilData = {};
    alltrialTimeStamps = {};
    alltrialSignalSummary = {};
    allMSData = {};
    counter = 1;
    for jj = 1 : 2 %2 %2 %2 % NOT 1:2 (for S07 / S08) 
        expcond = protocols{jj};
        sprintf(expcond)
        for kk = 1 : 8

            if ii==3 && jj == 1 && kk == 1
                continue
            elseif ii==7 && jj == 2 && kk == 2
                continue
            elseif ii==8 && jj==2 && kk==5 % this eyelink file is corrupted (missing data)
                continue
            end

            mdir = direction{kk};

            savePath = sprintf('%s/%s/ProcessedData/Summary/pupil', datadir, subj);

            if ~isfolder(savePath)
                mkdir(savePath)
            end
            
            MSpath = sprintf('%s/%s/ProcessedData/%s/eyedata/microsaccades/',datadir, subj, expcond);

            fileIDs = dir(sprintf('%s/%s/ProcessedData/%s/eyedata/MATs/%s*_tab.mat', datadir, subj, expcond, mdir));
            
            fileNames = {fileIDs.name}';
            fileNames = sort(fileNames);

            % look for MS event matrix
            mseventIDs = dir(sprintf('%s/%s/ProcessedData/%s/eyedata/microsaccades/*%s*_allevents.mat', datadir, subj, expcond, mdir));
            msfolder = mseventIDs(1).folder; % ok to index 1 - they are all in the same folder
            msfileNames = {mseventIDs.name}'; sort(msfileNames); msfileNames = flip(msfileNames); % this sorts with DIR2 first, so flip


            for ff=1:length(fileIDs)

                fileName = fileNames{ff}; % direction for that protocol

                load(sprintf('%s/%s/ProcessedData/%s/eyedata/MATs/%s', datadir, subj, expcond, fileName));
                sprintf('%s/%s/ProcessedData/%s/eyedata/MATs/%s', datadir, subj, expcond, fileName)
                
                if length(fileIDs)==1 || ff==1
                    pdata = load(sprintf('%s/%s/ProcessedData/%s/eyedata/components/%s%s_pupilMatrix.mat', datadir, subj, expcond, subj, mdir));
                    msdata = load(fullfile(msfolder, msfileNames{1}));
                    ms_path= fullfile(MSpath,sprintf('%s_microsaccadeMatrix.mat', direction{kk})); % 09/25
                else
                    % fix later: need to add the first block for those with 2 blocks 
                    pdata = load(sprintf('%s/%s/ProcessedData/%s/eyedata/components/%s%s2_pupilMatrix.mat', datadir, subj, expcond, subj, mdir));
                    msdata = load(fullfile(msfolder, msfileNames{2}));
                    ms_path= fullfile(MSpath,sprintf('%s2_microsaccadeMatrix.mat', direction{kk})); % new
                end

                try % in one case, eye tracker broke for first session, so only session 2 exists (S03) %new
                        load(ms_path);
                        ms_path
                catch
                        ms_path= fullfile(MSpath,sprintf('%s2_microsaccadeMatrix.mat', direction{kk}));
                        load(ms_path);
                        ms_path
                end

                % new MS log
                locations = 0:45:315; directions = 0:45:315; tilts = [0.5, 1, 2, 4, 8];
                [s_mean,RT_marker, RT_trialwise, RT_unfiltered, binaryRate, tabTrials] = computeMSRateRT(MS_TEMP,tab,1000, locations, directions, tilts, 4500);

                %msdata = msdata.eventMatrix==2; % 2 is MS label
                msdata = binaryRate;
                data = pdata.pupilMatrix;

                %msgLocate = dir(sprintf('%s/Data_DI_wEYE/%s/RawData/**/%s', datadir, subj, strrep(fileIDs.name, '_tab.mat', '.msg')));
                %msg_filepath = fullfile(msgLocate.folder, msgLocate.name);
                
                % correct RT for some early trials (did not account for stim period)
                RT_unfiltered = tab(:,13)*1000; % all RTs, no filter, convert to ms
                if any(RT_unfiltered<500) % if any trial RT is < 500 ms, that flags the early session needing the correction
                    RT_unfiltered = 500+RT_unfiltered;
                end
                
                [nTrials, nTimepoints] = size(data);
    
                %collectionRate = findSamplingRate(msg_filepath); % this is for anything using the timestamps (although everything else will be in 1000 hz)
                
                % change nans (at end of trial to 0)
                data(isnan(data)) = 0;
                
                % % initialize matrices
                % largetilt = [];
                % smalltilt = [];
                
                filteredPupilData = nan(size(data));
                trialTimeStamps = nan(nTrials,5);
                trialSignalSummary = nan(nTrials,4);
    
%                 if plotOn
%                     figure;
%                 end
                for di=1:nTrials
                    tempTrial = data(di,:);
                
                    lossTracking = find(abs(diff(tempTrial)) > (samplingRate/10) );
                    tempTrial(lossTracking) = 0;
                
                    % pad any 0s (no tracking / blinks) +- 10 ms with zeros
                    tempTrial = pad_zeros(tempTrial, pad_size);
                
                    tempTrial = blinkinterp(tempTrial,samplingRate,5,4,50,75);
                
                    tempTrial(isnan(tempTrial)) = 0;
                
                    if trialwise % define trialwise events
                
                        % retrieve trial-specific END
                        % I do not need to correct for sampling rate because this is in
                        % time units (e.g., 4555123, 4555125 (for a 500 hz - skips a value
                        % per ms)
                        trialEnd = (tab(di,8)-tab(di,2)); %/(collectionRate/samplingRate); % tab file time stamps needs to account for 1000 sampling
                        
                        % trial specific RT
                        trialResponse = RT_unfiltered(di)+1300; % here tab(:,13) is in seconds (RT) - no need to convert with collectionRate
                
                        % next trial Start
                        if di<nTrials
                            nextTrialStart = (tab(di+1,2)-tab(di,2)); %/(collectionRate/samplingRate); % tab file needs to account for 1000 sampling 
                        else
                            nextTrialStart = nan;
                        end
                
                
                        % remove end of trial 0s (for trialwise) - just cut the segment
                        % Find the index of the last nonzero element
                        last_nonzero_index = find(tempTrial ~= 0, 1, 'last');
                        
                        if isempty(last_nonzero_index)
                            % If the array is all zeros, no need to trim
                            trimmed_array = tempTrial;
                        else
                            % Trim the array to remove trailing zeros
                            trimmed_array = tempTrial(1:last_nonzero_index);
                        end
                
                        % if too many consecutive zeros, do not count in the fitting process 
                        % this means either eye was closed for too long, or next trial stimulus
                        % started too soon.
                
                    end
                
                    % smooth data
                    trimmed_array = myBWfilter(trimmed_array,[0.03 , 10],samplingRate,'bandpass');
                
                    if plotOn
                        figure(1)
                        plot(trimmed_array)
                        hold on
                        xline(stimStart, 'LineWidth',2)
                        hold on
                        xline(stimEnd, 'LineWidth',2)
                        hold on
                        xline(trialResponse, 'g')
                        hold on
                        xline(trialEnd, 'b')
                        hold on
                        if ~isnan(nextTrialStart)
                            xline(nextTrialStart, 'm')
                        end
                        hold off
                        xlim([1000 length(trimmed_array)]);
                        title(sprintf('Trial #%i', di))
                        hold on
                        %pause(.3)
                        %hold off
                    end
                
                    % all defined relative to trial start onset
                    trialTimeStamps(di,1) =stimStart; % in ms
                    trialTimeStamps(di,2) =stimEnd; % in ms
                    trialTimeStamps(di,3) =trialResponse; % in ms
                    trialTimeStamps(di,4) =trialEnd; % is ms
                    trialTimeStamps(di,5) =nextTrialStart; % in ms
                
                
                    filteredPupilData(di,1:length(trimmed_array)) = trimmed_array;
                



    %                 if tab(di,11)<4
    %                     %plot(tempTrial, 'Color',[.5 .5 .5])
    %                     smalltilt = [smalltilt; tempTrial];
    %                 else
    %                     %plot(tempTrial, 'k')
    %                     largetilt = [largetilt; tempTrial];
    %                 end
    %             %     hold on
    %             %     scatter(lossTracking, nanmean(tempTrial), 'r')
    %             %     title(sprintf('Trial #%i', di))
    %             %     hold on
    %             %     xline(1300)
    %             %     hold on
    %             %     xline(1800)
    %             %     hold on
    %             %     xline(trialEnd, 'r') % trial end % this needs to be adjusted based on samplingRate
    %             %     hold on
    %             %     xline(tab(di,13)*1000+1300, 'g') % this too (response time)
    %             %     hold on
    %             %     if di<nTrials
    %             %         xline(tab(di+1,2)-tab(di,2), 'r') 
    %             %     end
    %             %     
                
                end
    
                %  percent signal change
                
                tab = tab(1:nTrials,:); % use to determine if block ended early
                
                % Calculate the window mean/min and subtract that
                %meanPupilData = mean(filteredPupilData(:, 1000:1300), 2);
                %meanPupilData = meanPupilData(:);

                minPupilData = min(filteredPupilData(:, minWindowStart:minWindowEnd), [], 2); % note the [] to specify the dimension
                %minPupilData = min(filteredPupilData(:, 1000:1300), [], 2);
                [~, minIndices] = min(minPupilData, [], 2); % Find minimum value and index for each row
                minTimestamp = minWindowStart+minIndices-1; % log the trialwise min
                
                trialTimeStamps = [trialTimeStamps, minTimestamp]; % this is the minimum during stimulus period
                
                % THESE LINES BELOW WERE USED BEFORE
                % PSC: baseline does not take into account the stimulus
                meanPupilDataAllTime = abs(mean(filteredPupilData(:, baselinePSCstart:baselinePSCend), 2)); % to preserve the sign, take abs (within or across all sessions)?

                % first center the data so that prestimulus is always aligned (in size);
                % then divide by overall mean for that time period (within or across
                % sessions)? This converts to meaningful PSC
                %filteredPupilData = ((filteredPupilData - minPupilData)); % ./ nanmean(meanPupilDataAllTime) * 100; 
                

                % just trying
                %meanPupilDataAllTime = nanmean(filteredPupilData, "all");
                %filteredPupilData = (filteredPupilData - nanmin(nanmin(filteredPupilData))) ./ (nanmax(nanmax(filteredPupilData)) - nanmin(nanmin(filteredPupilData))) * 100; 
                %filteredPupilData = (filteredPupilData - meanPupilDataAllTime) ./ meanPupilDataAllTime * 100; 
                %minPupilData = min(filteredPupilData(:, minWindowStart:minWindowEnd), [], 2);
                minPupilData = mean(filteredPupilData(:, minWindowStart:minWindowEnd), 2);
                filteredPupilData = ((filteredPupilData - minPupilData));

                % just added 
                %filteredPupilData = (filteredPupilData - nanmean(nanmean(filteredPupilData))) ./ nanmean(nanmean(filteredPupilData)) * 100;
                
                %disp('MEAN')
                %nanmean(meanPupilDataAllTime)
                % THESE LINES ABOVE WERE USED BEFORE

% NEW WAY

                [nTrials,~] = size(filteredPupilData);
                aa = filteredPupilData; %(:, 1000:1800); %1300);
                amaxtemp = max(aa'); amax = max(amaxtemp);
                amintemp = min(aa'); amin = min(amintemp);
                %rangePupilDataAllTime = max(abs(amax-amin)); % range across trials
                rangePupilDataAllTime = amaxtemp-amintemp; 
                %rangePupilDataAllTime = rangePupilDataAllTime';
                rangePupilDataAllTime = repmat(mean(rangePupilDataAllTime), nTrials, 1); 

%                 % new way (compute range - 
%                 aa = filteredPupilData; %(:, 1000:1800); %1300);
%                 amaxtemp = max(aa'); amax = max(amaxtemp);
%                 amintemp = min(aa'); amin = min(amintemp);
%                 rangePupilDataAllTime = max(abs(amax-amin)); % range across trials (per session)

%                 aa = filteredPupilData(:, 1000:1800); %1300);
%                 amaxtemp = mean(aa'); amax = max(amaxtemp);
%                 amintemp = mean(aa'); amin = min(amintemp);
%                 rangePupilDataAllTime = amax-amin;
%                 %rangePupilDataAllTime = amaxtemp-amintemp; %rangePupilDataAllTime = rangePupilDataAllTime'; % trialwise range
%                 %disp('Range')
%                 %rangePupilDataAllTime
% 
%                 %maxPupilData = max(filteredPupilData(:, minWindowStart:minWindowEnd), [], 2); % note the [] to specify the dimension
%                 %[~, maxIndices] = max(minPupilData, [], 2); % Find minimum value and index for each row
%                 %max(maxPupilData-minPupilData)
%                 signalDiff = filteredPupilData - minPupilData;
%                 filteredPupilData = ((signalDiff ./ rangePupilDataAllTime) * 100);
% 
%                 
%                 %meanAllSignal = abs(nanmean(signalDiff, 'all'));
%                 %filteredPupilData = (signalDiff ./ meanAllSignal) * 100; 

                trialSignalSummary(:,1) = minPupilData;
                %trialSignalSummary(:,2) = meanPupilDataAllTime;
                trialSignalSummary(:,2) = rangePupilDataAllTime;
                trialSignalSummary(:,3) = mean(meanPupilDataAllTime); %nanmax(nanmax(filteredPupilData)); %nanmedian(nanmedian(filteredPupilData)); %
                %trialSignalSummary(:,3) = mean(rangePupilDataAllTime);
                trialSignalSummary(:,4) = meanPupilDataAllTime;

%                 % remove later
%                 if plotOn
%                     for di=1:nTrials
%                         plot(filteredPupilData(di,:))
%                         hold on
%                         xline(stimStart, 'LineWidth',2)
%                         hold on
%                         xline(stimEnd, 'LineWidth',2)
%                         hold on
%                         xline(trialResponse, 'g')
%                         hold on
%                         xline(trialEnd, 'b')
%                         hold on
%                         if ~isnan(nextTrialStart)
%                             xline(nextTrialStart, 'm')
%                         end
%                         hold off
%                         xlim([1000 length(trimmed_array)]);
%                         title(sprintf('Trial #%i', di))
%                         pause(.3)
%                         hold off
%                     end
%                 end



    
%                 if plotOn
%                     figure(2)
%                     hardTrials = filteredPupilData(tab(:,11)<=4,:);
%                     easyTrials = filteredPupilData(tab(:,11)>4,:);
%                     
%                     RT_hard = trialTimeStamps(tab(:,11)<=4,3); 
%                     RT_easy = trialTimeStamps(tab(:,11)>4,3); 
%                     
%                     medianRT_hard = mean(RT_hard(RT_hard<2000));
%                     medianRT_easy = mean(RT_easy(RT_easy<2000));
%                     
%                     plot(nanmedian(hardTrials,1), 'Color', [.5 .5 .5])
%                     hold on
%                     plot(nanmedian(easyTrials,1), 'k')
%                     hold on
%                     xline(1300, 'LineWidth',2)
%                     hold on
%                     xline(1800, 'LineWidth',2)
%                     hold on
%                     xline(medianRT_hard, 'r')
%                     hold on
%                     xline(medianRT_easy, 'g')
% 
% 
%                     % card / oblique
%                     figure(3)
%                     easyTrialsIdx = find(ismember(tab(:,10), 5:8)); % & tabnew(:,14)==1; 
%                     hardTrialsIdx = find(ismember(tab(:,10), 1:4)); % & tabnew(:,14)==1;
%                     
%                     hardTrials = filteredPupilData(hardTrialsIdx,:);
%                     easyTrials = filteredPupilData(easyTrialsIdx,:);
%                     
%                     RT_hard = trialTimeStamps(hardTrialsIdx,3); 
%                     RT_easy = trialTimeStamps(easyTrialsIdx,3); 
%                     
%                     medianRT_hard = mean(RT_hard(RT_hard<2000));
%                     medianRT_easy = mean(RT_easy(RT_easy<2000));
%                     
%                     plot(nanmedian(hardTrials,1), 'Color', [51, 34, 136]/255)
%                     hold on
%                     plot(nanmedian(easyTrials,1), 'Color', [17, 119, 51]/255)
%                     hold on
%                     xline(1300, 'LineWidth',2)
%                     hold on
%                     xline(1800, 'LineWidth',2)
%                     hold on
%                     xline(medianRT_hard, 'r')
%                     hold on
%                     xline(medianRT_easy, 'g')
% 
%                 end
    
                allfilteredPupilData{counter} = filteredPupilData;
                alltrialTimeStamps{counter} = trialTimeStamps;
                alltab{counter} = tab;
                alltrialSignalSummary{counter} = trialSignalSummary;
                allMSData{counter} = msdata;

                counter = counter+1;
            end
        end
    end

    % Unravel data from all sessions
    allfilteredPupilData = vertcat(allfilteredPupilData{:});
    alltrialTimeStamps = vertcat(alltrialTimeStamps{:});
    alltab = vertcat(alltab{:});
    alltrialSignalSummary = vertcat(alltrialSignalSummary{:});

    % special case for S03- data truncated for cell 11 and 14
    if ii == 3
        allMSData{11} = allMSData{11}(1:521,:);
        allMSData{14} = allMSData{14}(1:644,:);
    end

    allMSData = vertcat(allMSData{:});

    % put back!!
    %allfilteredPupilData = (allfilteredPupilData ./ max(alltrialSignalSummary(:,2))) * 100; % divide by max range
    %allfilteredPupilData = (allfilteredPupilData ./ alltrialSignalSummary(:,2)) * 100; 

    %deno = (max(alltrialSignalSummary(:,2))*(max(alltrialSignalSummary(:,2))./(alltrialSignalSummary(:,2))));
    %allfilteredPupilData = (allfilteredPupilData ./ alltrialSignalSummary(:,2)) * 100; % divide by max range

    allfilteredPupilData = (allfilteredPupilData - nanmean(nanmean(allfilteredPupilData))) ./ abs(nanmean(nanmean(allfilteredPupilData)));

    if plotOn
        figure
        hardTrials = allfilteredPupilData(alltab(:,11)<=4,:);
        easyTrials = allfilteredPupilData(alltab(:,11)>4,:);
        
        RT_hard = alltrialTimeStamps(alltab(:,11)<=4,3); 
        RT_easy = alltrialTimeStamps(alltab(:,11)>4,3); 
        
        medianRT_hard = mean(RT_hard(RT_hard<2000));
        medianRT_easy = mean(RT_easy(RT_easy<2000));
        
        plot(nanmean(hardTrials,1), 'Color', [.5 .5 .5])
        hold on
        plot(nanmedian(easyTrials,1), 'k')
        hold on
        xline(1300, 'LineWidth',2)
        hold on
        xline(1800, 'LineWidth',2)
        hold on
        xline(medianRT_hard, 'r')
        hold on
        xline(medianRT_easy, 'g')


        % card / oblique
        figure
        easyTrialsIdx = find(ismember(alltab(:,10), 5:8)); % & tabnew(:,14)==1; 
        hardTrialsIdx = find(ismember(alltab(:,10), 1:4)); % & tabnew(:,14)==1;
        
        hardTrials = allfilteredPupilData(hardTrialsIdx,:);
        easyTrials = allfilteredPupilData(easyTrialsIdx,:);
        
        RT_hard = alltrialTimeStamps(hardTrialsIdx,3); 
        RT_easy = alltrialTimeStamps(easyTrialsIdx,3); 
        
        medianRT_hard = mean(RT_hard(RT_hard<2000));
        medianRT_easy = mean(RT_easy(RT_easy<2000));
        
        plot(nanmedian(hardTrials,1), 'Color', [51, 34, 136]/255)
        hold on
        plot(nanmedian(easyTrials,1), 'Color', [17, 119, 51]/255)
        hold on
        xline(1300, 'LineWidth',2)
        hold on
        xline(1800, 'LineWidth',2)
        hold on
        xline(medianRT_hard, 'r')
        hold on
        xline(medianRT_easy, 'g')

        figure
        subplot(1,3,1)
        scatter(ones(length(easyTrialsIdx),1), alltrialSignalSummary(easyTrialsIdx,2), 'b')
        hold on
        scatter(1, mean(alltrialSignalSummary(easyTrialsIdx,2)),300, 'b+', 'LineWidth', 3)
        hold on
        scatter(1+ones(length(hardTrialsIdx),1), alltrialSignalSummary(hardTrialsIdx,2), 'r')
        hold on
        scatter(2, mean(alltrialSignalSummary(hardTrialsIdx,2)),300, 'r+', 'LineWidth', 3)
        title('range')
        xlim([0, 3])
        subplot(1,3,2)
        scatter(ones(length(easyTrialsIdx),1), alltrialSignalSummary(easyTrialsIdx,3), 'b')
        hold on
        scatter(1, mean(alltrialSignalSummary(easyTrialsIdx,3)),300, 'b+', 'LineWidth', 3)
        hold on
        scatter(1+ones(length(hardTrialsIdx),1), alltrialSignalSummary(hardTrialsIdx,3), 'r')
        hold on
        scatter(2, mean(alltrialSignalSummary(hardTrialsIdx,3)),300, 'r+', 'LineWidth', 3)
        xlim([0, 3])
        title('session mean')
        subplot(1,3,3)
        scatter(ones(length(easyTrialsIdx),1), alltrialSignalSummary(easyTrialsIdx,4), 'b')
        hold on
        scatter(1, mean(alltrialSignalSummary(easyTrialsIdx,4)),300, 'b+', 'LineWidth', 3)
        hold on
        scatter(1+ones(length(hardTrialsIdx),1), alltrialSignalSummary(hardTrialsIdx,4), 'r')
        hold on
        scatter(2, mean(alltrialSignalSummary(hardTrialsIdx,4)),300, 'r+', 'LineWidth', 3)
        title('trial mean')
        xlim([0, 3])

    end


    % save it out
    save(fullfile(savePath,sprintf('%s_allpupilData.mat', subj)), 'allfilteredPupilData', 'alltrialTimeStamps', 'alltab', 'alltrialSignalSummary');
    save(fullfile(strrep(savePath, 'pupil', 'microsaccades'), sprintf('%s_allmsData.mat', subj)), 'allMSData');
end

%% model fitting

% for modeling
samplingRate = 1000; % hertz
cutOff = 1300; % this is to make the time meaningful (STIM ON = 0)
event_label = {'pre_stim', 'sti on', 'response'}; % events modeled (other than boxcar)
condlabels = {'trial'}; % arbitrary
baseline = []; % empty because I did the normalization myself a different way    
wnum = 1;

for ii =4:4
    subj = subjects{ii};
    savePath = sprintf('%s/%s/ProcessedData/Summary/pupil', datadir, subj);
    load(fullfile(savePath,sprintf('%s_allpupilData.mat', subj)));

    [nTrials, ~] = size(allfilteredPupilData);
    
    field_names = {'eventtimes', 'boxtimes', 'samplerate', 'window', 'ampvals', 'boxampvals', 'latvals', 'tmaxval', 'yintval', 'slopeval','numparams','cost','R2','BICrel'};  % Add more field names as needed
    
    % Preallocate the struct array
    for i = 1:numel(field_names)
        [output(1:nTrials) .(field_names{i})] = deal([]);  % Initialize each field with empty values
    end
    
    tic
    
    parfor trialNum=1:nTrials
    %for trialNum=35 %20 %12 %1:nTrials
    
        if mod(trialNum, 100) == 0
            sprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Trial# %i ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', trialNum)
        end
    
        tempData = allfilteredPupilData(trialNum,:); %nanmedian(hardTrials,1);
        
    %     if tab(trialNum,11)<4
    %         disp('HARD TRIAL')
    %     else
    %         disp('EASY TRIAL')
    %     end
        
        timestamps = alltrialTimeStamps(trialNum,:);
        
        TrialEnd = timestamps(4); %alltrialTimeStamps(trialNum,4); %nanmedian(trialTimeStamps(:,4), 1);
        responseTime = timestamps(3); %alltrialTimeStamps(trialNum,3);
        minDuringStim = timestamps(6); %alltrialTimeStamps(trialNum,6);
        
        % fit the model to the mean for a subject
        % figure
        % plot(tempData)
        
        % event_label = {'fix enforced', 'sti on', 'response cue', 'response', 'next trial start'};
        % event = [1000-cutOff, 1300-cutOff, 1800-cutOff, responseTime-cutOff, medianTrialEnd-cutOff];
        
        event = [800-cutOff, stimStart-cutOff, responseTime-cutOff];
        
        trialwindow = [-cutOff length(tempData)-(cutOff+1)]; % define the timestamps? 0 = stim Onset
        %modelwindow = [cutOff length(tempData)];
    
        options = pret_preprocess();
        options.normflag = false;
        options.blinkflag = false;
        
        sj = pret_preprocess({tempData},samplingRate,trialwindow,condlabels,baseline,options);
        
        % create model for the data
        % Pretending that we are naive to the model that we used to create our
        % data, let's create a model to actually fit the data to.
        model = pret_model();
        
        % While the trial window of our task is from -500 to 3500 ms, here we are
        % not interested in what's happening before 0. So
        % let's set the model window to fit only to the region betweeen 0 and 3500
        % ms (the cost function will only be evaluated along this interval).
        model.window = [800-cutOff trialwindow(2)]; % minDuringStim-cutOff        %[-800 trialwindow(2)];
        
        % We already know the sampling frequency.
        model.samplerate = samplingRate;
        
        % We also know the event times of our task. Let's also say that we think 
        % there will be a sustained internal signal from precue onset to response 
        % time (0 to 2750 ms).
        model.eventtimes = event;
        model.eventlabels = event_label; %optional
        model.boxtimes = {[stimStart-cutOff stimEnd-cutOff]}; %, [-300 500]}; %round(responseTime-1300)]}; %round(responseTime-1300)]}; % should I make boxcar task or stimulus related?
        model.boxlabels = {'stim'}; %, 'fixation'}; %optional
        
        % Let's say we want to fit a model with the following parameters: 
        % event-related, amplitude, latency, task-related (box) amplitude, 
        % and the tmax of the pupil response function. We turn the other parameters
        % off.
        model.yintflag = false;
        model.slopeflag = false;
        model.tmaxflag = false;
        % Now let's define the bounds for the parameters we decided to fit. We do
        % not have to give values for the y-intercept and slope because we are not
        % fitting them.
        model.ampbounds = repmat([0;200],1,length(model.eventtimes)); % changed to 150 from 100 b/c optimization frequency hits the max
        %model.latbounds = repmat([-500;500],1,length(model.eventtimes));
        %model.latbounds = [-1000 0 0; 0 responseTime-stimStart, 1000];
        model.latbounds = [-1000 0 -500; 500 responseTime-stimStart, 1500];

        %event_label = {'fix enforced', 'sti on', 'response cue', 'response', 'next trial start'};
        model.boxampbounds = [0;200]; %[0 0 ;100 100];
        %model.tmaxbounds = [800;800];
        
        % We need to fill in the values for the y-intercept and slope since we will
        % not be fitting them as parameters.
        model.yintval = 0;
        model.slopeval = 0;
        model.tmaxval = 930; %800;
        
        % estimate model parameters via pret_estimate_sj
        % Now let's perform the parameter estimation procedure on our subject data.
        % The mean of each condition will be fit independently. For illustration, 
        % let's run only 3 optimizations using one cpu worker (for more 
        % information, see the help files of pret_estimate and pret_estimate_sj).
        options = pret_estimate_sj();
        options.pret_estimate.optimnum = 3;
        % if you want to try fiting the parameters using single trials instead of the mean,
        % use these lines (you'll want to turn off the optimization plots for this):
        %options.trialmode = 'single';
        
        options.pret_estimate.pret_optim.optimplotflag = false; %true; %
        
        sj = pret_estimate_sj(sj,model,wnum,options);
    
    %     % to add peak to sj - no longer needed: always 930
    %     % Y calc is the sum, X: row per event - each representing a gaussian
    %     sfact = samplingRate/1000;
    %     time = model.window(1):1/sfact:model.window(2);
    %     [Ycalc, X] = pret_calc(sj.estim.trial); % use this to extract peak
    %     X = X(1:end-1,:);
    %     [max_values, max_indices] = max(X, [], 2); 
    % 
    %     % convert to time (relative to latency)
    %     peakmsfromZero = time(max_indices);
    %     peakmsfromLatency = peakmsfromZero - (sj.estim.trial.eventtimes+sj.estim.trial.latvals); % 930
    
        %output = [output; sj.estim.trial];
    
        output(trialNum) = sj.estim.trial;
    
    end
    toc
    
    % save it out
    save(fullfile(savePath,sprintf('%s_allpupilFits.mat', subj)), 'output');
    delete(gcp); % shut down parallel pool
end

%% helpers

function padded_array = pad_zeros(array, pad_size)
    padded_array = array; % Make a copy of the original array
    zero_indices = find(array == 0); % Find indices of zeros
    
    % Iterate over zero indices
    for idx = zero_indices
        start_pad = max(1, idx - pad_size); % Calculate start padding index
        end_pad = min(length(array), idx + pad_size); % Calculate end padding index
        
        % Pad with zeros if padding is within array boundaries
        padded_array(start_pad:idx-1) = 0;
        padded_array(idx+1:end_pad) = 0;
    end
end

