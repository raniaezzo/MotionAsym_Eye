clc; clear all; %close all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
%datadir='/Users/rje257/Desktop/MS_Project/Data_DI_wEYE'

analysis_type = 'direction';

stimonset = 1300; % this data has 1000 ms cutoff from beginning
stimoffset = 1800; 

analyses = {'direction', 'tilt', 'outcome'};
x_values = linspace(-5, 5, 100);
z_score_init = 60000; % arbitrary, but never goes beyond

msOnset = nan(length(subjects), 2, length(analyses)); % subjects, easy/hard, analysis
zscored_msOnset = nan(length(subjects), 2, length(analyses));
plotOn = 0;

% Initialize the cell array for probability density
probDensity = cell(1, length(analyses));
grand_z_scores_MSonset = cell(1, length(analyses)); 
grand_possible = cell(1, length(analyses)); 

amp = cell(1, length(analyses));
boxamp = cell(1, length(analyses));
latvals = cell(1, length(analyses));
eventtimes = cell(1, length(analyses));

for cc = 1:length(analyses)
    matrix = nan(length(subjects), length(x_values), 2);
    matrix2 = nan(length(subjects), z_score_init, 2);
    matrix3 = nan(length(subjects), 1, 2);
    matrix4 = nan(length(subjects),3,2);
    % Assign the matrix to the cell
    probDensity{cc} = matrix;
    grand_z_scores_MSonset{cc} = matrix2; % mat2 is for zscores
    grand_possible{cc} = matrix3; % mat2 is for zscores
    amp{cc} = matrix4;
    boxamp{cc} = matrix3; % only 1 val
    latvals{cc} = matrix4;
    eventtimes{cc} = matrix4;
end

totalQualifiedTrials = nan(length(subjects), length(analyses));

grand_z_scores_MSonset_easy = [];
grand_z_scores_MSonset_hard = [];
grand_possible_easy = [];
grand_possible_hard = [];

rtConds = nan(length(subjects), 3, length(analyses)); % 3 time measures saved

if plotOn ~= 0
    %f1 = figure(1);
    f2 = figure(2);
    f3 = figure(3);
    f4 = figure(4);
    f5 = figure(5);
    f6 = figure(6);
end

f6.Position = [135 739 1848 400];

% get the scale factor for PSC
psc_scaleFactor = nan(length(subjects), 1);
    
%%
for si=1:length(subjects)

    f1 = figure;
    subj=subjects{si};

    savePath = sprintf('%s/%s/ProcessedData/Summary/pupil', datadir, subj);
    load(fullfile(savePath,sprintf('%s_allpupilFits.mat', subj)), 'output');
    load(fullfile(savePath,sprintf('%s_allpupilData.mat', subj))); % to load in tab file

    psc_scaleFactor(si,1) = (nanmean(alltrialSignalSummary(:,2))/nanmax(alltrialSignalSummary(:,4)));


    rate = load(fullfile(strrep(savePath, 'pupil', 'microsaccades'), sprintf('%s_allmsData.mat', subj)));

    % S07 has VL file missing
    if si==7
         alltrialTimeStamps = alltrialTimeStamps(alltab(:,10)~=7,:);
         alltrialSignalSummary = alltrialSignalSummary(alltab(:,10)~=7,:);
         allfilteredPupilData = allfilteredPupilData(alltab(:,10)~=7,:);
         alltab = alltab(alltab(:,10)~=7,:);
         rate.allMSData = rate.allMSData(alltab(:,10)~=7,:);
    end

%     % get rows that are cardinal / oblique 
%     easyTrials = ismember(alltab(:,10), 5:8); 
%     hardTrials = ismember(alltab(:,10), 1:4);
% 
%     rate_easy = sum(rate.allMSData(easyTrials,:),1);
%     rate_hard = sum(rate.allMSData(hardTrials,:),1);

    % find the mean suppression point (during stimulus)
    %combinedData = nanmean([rate.(fieldNames{1}); rate.(fieldNames{2})]);
    %[~,suppressionTime] = min(combinedData(300:800)); % stim on to offset (adds back 300 in the function)

%     meanMS = nanmean(rate.allMSData,1);
%     [~,suppressionTime] = min(meanMS(stimonset:stimoffset)); 

    % Correction for some trials with <500 RT: some early sessions counted
    % ms relative to stimOff rather than stimOn
    temp = alltab(:,1); % get the trial number
    next_session = [false; diff(temp) < 0];     % Find the points where the numbers start to decrease
    indices = cumsum(next_session) + 1;
    rtTrialwise.allMSData = nan(length(alltab),1);
    suppressionTimes = nan(length(alltab),1);
    postsupPeakTimes = nan(length(alltab),1);

    % to find the post stim inflection point
    window_size = 250; % Window size for the moving average

    for sn=1:max(indices) % this is by session
        [selectTrials, ~] = find(indices==sn);
    
%         % find the suppression for the mean rate over session
        meanMS = nanmean(rate.allMSData(selectTrials,:),1);

        smoothed_data = smoothdata(meanMS, 'movmean', window_size);
            
        decreasing_index = find(diff(smoothed_data(1800:end)) < 0, 1);

        if isempty(decreasing_index)
            decreasing_index = find(diff(smoothed_data(1800:end)) == 0, 1); % if the MS rate does not go back down, select when it plateaus
        end

        postsupPeakTimes(selectTrials) = ones(length(selectTrials),1)*(decreasing_index+1800);

        [~,suppressionTime] = min(meanMS(stimonset:stimoffset)); 
        suppressionTimes(selectTrials) = ones(length(selectTrials),1)*suppressionTime; % suppression computed within session (but this doesn't matter)

        RT_unfiltered = alltab(selectTrials,13)*1000; % all RTs, no filter
        if any(RT_unfiltered<500)
            rtTrialwise.allMSData(selectTrials) = RT_unfiltered +500;
        else
            rtTrialwise.allMSData(selectTrials) = RT_unfiltered;
        end
    end

    rtTrialwise.allMSData(rtTrialwise.allMSData>=2000) = nan;

    trialData.allMSData = alltab;

    [postOnsetMS_cond, latencyCond, RT_outputCond, accuracy_Cond, qualifyingTrials, dataCheck, supTimes, postOffsetTimes, selectedTrials] = ...
        compute_latency(rate, 'allMSData', suppressionTimes, postsupPeakTimes, rtTrialwise, trialData, stimonset);

    postOnsetMS_cond = postOnsetMS_cond; % - 1300; % to give 0 stimonset meaning
    %postOnsetMS_cond = latencyCond;

    tabnew = qualifyingTrials{1}; % 1 b/c all conditions are in cell 1

    pupilFitsMatch = output(selectedTrials);
    pupilData = allfilteredPupilData(selectedTrials,:);
    
    for ai=1:length(analyses)

        analysis_type = analyses{ai};

        if plotOn ~= 0
            figure(1)
            subplot(1,3,ai)
        end

        % Correalte MS latency with Pupil Latency
        if strcmp(analysis_type, 'direction')
            fieldNames = {'cardinal', 'oblique'};
            color = {[17, 119, 51],[51, 34, 136]}; 
        elseif strcmp(analysis_type, 'tilt')
            fieldNames = {'largeoffset','smalloffset'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'location')
            fieldNames = {'horizontalLoc','verticalLoc'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'dirtilt')
            fieldNames = {'easycardinal','hardoblique'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'outcome')
            fieldNames = {'correct','incorrect'};
            color = {[0, 255, 0],[255, 0, 0]}; 
        end


        if strcmp(analysis_type, 'direction')
            easyTrials = ismember(tabnew(:,10), 5:8); % & tabnew(:,14)==1; 
            hardTrials = ismember(tabnew(:,10), 1:4); % & tabnew(:,14)==1;
            possibleEasy = sum(ismember(trialData.allMSData(:,10), 5:8));
            possibleHard = sum(ismember(trialData.allMSData(:,10), 1:4));
        elseif strcmp(analysis_type, 'tilt')
%             easyTrials = tabnew(:,11)>=4; %  & tabnew(:,14)==1;
%             hardTrials = tabnew(:,11)<2; %  & tabnew(:,14)==1;
%             possibleEasy = sum(trialData.allMSData(:,11)>4); 
%             possibleHard = sum(trialData.allMSData(:,11)<=4); %2);
            easyTrials = tabnew(:,11)==8; %  & tabnew(:,14)==1;
            hardTrials = tabnew(:,11)<=4; %  & tabnew(:,14)==1;
            possibleEasy = sum(trialData.allMSData(:,11)==8); 
            possibleHard = sum(trialData.allMSData(:,11)<=4); %2);
        elseif strcmp(analysis_type, 'outcome')
            easyTrials = tabnew(:,14)==1;
            hardTrials = tabnew(:,14)==0;
            possibleEasy = sum(trialData.allMSData(:,14)==1);
            possibleHard = sum(trialData.allMSData(:,14)==0);
        end

        if plotOn ~= 0
            scatter(RT_outputCond(easyTrials), postOnsetMS_cond(easyTrials), 'o', 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
            hold on
            scatter(RT_outputCond(hardTrials), postOnsetMS_cond(hardTrials), 'o', 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
            hold on
            scatter(nanmedian(RT_outputCond(easyTrials)), nanmean(postOnsetMS_cond(easyTrials)), 350, 'o', 'filled', 'MarkerFaceColor', [1 1 1])
            hold on
            scatter(nanmedian(RT_outputCond(easyTrials)), nanmean(postOnsetMS_cond(easyTrials)), 350, '+', 'MarkerEdgeColor', color{1}/255, 'LineWidth', 5, 'MarkerFaceAlpha', 0)
            hold on
            scatter(nanmedian(RT_outputCond(hardTrials)), nanmean(postOnsetMS_cond(hardTrials)), 350, 'o', 'filled', 'MarkerFaceColor', [1 1 1])
            hold on
            scatter(nanmedian(RT_outputCond(hardTrials)), nanmean(postOnsetMS_cond(hardTrials)), 350, '+', 'MarkerEdgeColor', color{2}/255, 'LineWidth', 5, 'MarkerFaceAlpha', 0)
            xlim([500 2000])
            xlabel('reaction time (ms)')
            ylabel('onset of first microsaccade (ms)')
            title(analysis_type)
            %ylabel('latency (ms)')
        end

        % correlation_coefficient1 = corr(RT_outputCond1, latencyCond1', 'type', 'Spearman')
        correlation_coefficient1 = corr(RT_outputCond(easyTrials), postOnsetMS_cond(easyTrials)', 'type', 'Spearman', 'rows', 'complete')
        correlation_coefficient2 = corr(RT_outputCond(hardTrials), postOnsetMS_cond(hardTrials)', 'type', 'Spearman', 'rows', 'complete')
    
        msOnset(si,1, ai) = nanmean(postOnsetMS_cond(easyTrials));
        msOnset(si,2, ai) = nanmean(postOnsetMS_cond(hardTrials));

        if plotOn ~= 0
            [x_fit, y_fit] = fitLine(RT_outputCond(easyTrials), postOnsetMS_cond(easyTrials));
            plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', color{1}/255);
        
            hold on
            [x_fit2, y_fit2] = fitLine(RT_outputCond(hardTrials), postOnsetMS_cond(hardTrials));
            plot(x_fit2, y_fit2, '-', 'LineWidth', 2, 'Color', color{2}/255);
    
            text(1800, 1200, num2str(correlation_coefficient1)); % Coordinates for the first line of text (2, 5)
            text(1800, 1100, num2str(correlation_coefficient2)); 
        end

        % recalculate the suppression for the hard vs. easy (based on data output)
        
        [postsupPeakTimes_easy, supTime_easy] = calculateSuppPeak(rate.allMSData, easyTrials, window_size);
        [postsupPeakTimes_hard, supTime_hard] = calculateSuppPeak(rate.allMSData, hardTrials, window_size);
        
        rtConds(si,1,ai) = nanmedian(RT_outputCond(easyTrials));
        rtConds(si,2,ai) = nanmedian(RT_outputCond(hardTrials));
        rtConds(si,3,ai) = nanmedian(RT_outputCond(hardTrials))-nanmedian(RT_outputCond(easyTrials));

        if plotOn ~= 0
            figure(2)
            subplot(1,3,ai)
            %scatter(ones(length(postOffsetTimes(easyTrials)),1), postOffsetTimes(easyTrials), 'filled', 'MarkerFaceColor', color{1}/255)
            scatter(1, supTime_easy, 'filled', 'MarkerFaceColor', color{1}/255)
            hold on
            %scatter(2*ones(length(postOffsetTimes(hardTrials)),1), postOffsetTimes(hardTrials), 'filled', 'MarkerFaceColor', color{2}/255)
            scatter(2, supTime_hard, 'filled', 'MarkerFaceColor', color{2}/255)
            hold on
            plot([1,2], [supTime_easy, supTime_hard], 'k')
            title('Suppression')
            xlim([0 3])
    
            figure(3)
            subplot(1,3,ai)
            %scatter(nanmedian(RT_outputCond(easyTrials)), supTime_easy, 'filled', 'MarkerFaceColor', color{1}/255)
            scatter(1, nanmedian(RT_outputCond(easyTrials)), 'filled', 'MarkerFaceColor', color{1}/255)
            hold on
            scatter(2, nanmedian(RT_outputCond(hardTrials)), 'filled', 'MarkerFaceColor', color{2}/255)
            hold on
            plot([1,2], [nanmedian(RT_outputCond(easyTrials)), nanmedian(RT_outputCond(hardTrials))])
            %scatter(2*ones(length(postOffsetTimes(hardTrials)),1), postOffsetTimes(hardTrials), 'filled', 'MarkerFaceColor', color{2}/255)
            %scatter(nanmedian(RT_outputCond(hardTrials)), supTime_hard, 'filled', 'MarkerFaceColor', color{2}/255)
            title('RT')
            xlim([0 3])
    
    %         scatter(1, postsupPeakTimes_easy, 'filled', 'MarkerFaceColor', color{1}/255)
    %         hold on
    %         %scatter(2*ones(length(postOffsetTimes(hardTrials)),1), postOffsetTimes(hardTrials), 'filled', 'MarkerFaceColor', color{2}/255)
    %         scatter(2, postsupPeakTimes_hard, 'filled', 'MarkerFaceColor', color{2}/255)
    %         hold on
    %         plot([1,2], [postsupPeakTimes_easy, postsupPeakTimes_hard], 'k')
    %         title('Post Suppression Peak')
    %         xlim([0 3])
    
            figure(4)
            subplot(1,3,ai)
            scatter(1, postsupPeakTimes_easy-supTime_easy, 'filled', 'MarkerFaceColor', color{1}/255)
            hold on
            %scatter(2*ones(length(postOffsetTimes(hardTrials)),1), postOffsetTimes(hardTrials), 'filled', 'MarkerFaceColor', color{2}/255)
            scatter(2, postsupPeakTimes_hard-supTime_hard, 'filled', 'MarkerFaceColor', color{2}/255)
            hold on
            plot([1,2], [postsupPeakTimes_easy-supTime_easy, postsupPeakTimes_hard-supTime_hard], 'k')
            title('Duration of Min to Max')
            xlim([0 3])
    
    
            figure(5)
            subplot(1,3,ai)
            scatter(nanmedian(RT_outputCond(hardTrials))-nanmedian(RT_outputCond(easyTrials)), supTime_hard-supTime_easy, 'filled', 'MarkerFaceColor', 'k')
            nanmedian(RT_outputCond(hardTrials))-nanmedian(RT_outputCond(easyTrials))
            supTime_hard-supTime_easy
            hold on
            title('RT diff vs Suppr diff')
        end

        if plotOn == 0
            figure(6)
            subplot(1,3,ai)
            plot(nanmean(pupilData(easyTrials, :),1), 'Color', color{1}/255, 'LineWidth', 2)
            hold on
            plot(nanmean(pupilData(hardTrials, :),1), 'Color', color{2}/255, 'LineWidth', 2)
            hold on
            xlim([1300 4500])
        end

        % pupil data
        pupilFit_easy = pupilFitsMatch(easyTrials);
        pupilFit_hard = pupilFitsMatch(hardTrials);

        varExp_easy = cell2mat({pupilFit_easy.('R2')});
        varExp_hard = cell2mat({pupilFit_hard.('R2')});

        % which values do I need to plot the signal (reconstruction)?
        pupilFit_easy_filt = pupilFit_easy(varExp_easy>0.50);
        pupilFit_hard_filt = pupilFit_hard(varExp_hard>0.50);

        amp{ai}(si,1:3,1) = median(cell2mat({pupilFit_easy_filt.('ampvals')}'), 1);
        boxamp{ai}(si,1,1) = median(cell2mat({pupilFit_easy_filt.('boxampvals')}'), 1);
        latvals{ai}(si,1:3,1) = median(cell2mat({pupilFit_easy_filt.('latvals')}'), 1);
        eventtimes{ai}(si,1:3,1) = median(cell2mat({pupilFit_easy_filt.('eventtimes')}'), 1);
        amp{ai}(si,1:3,2) = median(cell2mat({pupilFit_hard_filt.('ampvals')}'), 1);
        boxamp{ai}(si,1,2) = median(cell2mat({pupilFit_hard_filt.('boxampvals')}'), 1);
        latvals{ai}(si,1:3,2) = median(cell2mat({pupilFit_hard_filt.('latvals')}'), 1);
        eventtimes{ai}(si,1:3,2) = median(cell2mat({pupilFit_hard_filt.('eventtimes')}'), 1);



        % take the z-score to allow for combinging data across subjects
        %z_scores_RT = zscore(RT_outputCond);

        avTotal = mean(postOnsetMS_cond);
        stdevTotal = std(postOnsetMS_cond);
        %z_scored_data = (postOnsetMS_cond - avTotal) / stdevTotal;
        z_scored_data = (postOnsetMS_cond-1300); % / stdevTotal) * -1; % multiply by -1 to make negative interpretation as "early time"

        z_scores_MSonset_easy = z_scored_data(easyTrials);
        z_scores_MSonset_hard = z_scored_data(hardTrials);

        grand_z_scores_MSonset{ai}(si,1:length(z_scores_MSonset_easy),1) = z_scores_MSonset_easy;
        grand_z_scores_MSonset{ai}(si,1:length(z_scores_MSonset_hard),2) =  z_scores_MSonset_hard;
        grand_possible{ai}(si,1:length(possibleEasy),1) =  possibleEasy;
        grand_possible{ai}(si,1:length(possibleHard),2) =  possibleHard;

%         z_scores_MSonset_easy = zscore(postOnsetMS_cond(easyTrials));
%         z_scores_MSonset_hard = zscore(postOnsetMS_cond(hardTrials));

        zscored_msOnset(si,1, ai) = nanmedian(z_scores_MSonset_easy);
        zscored_msOnset(si,2, ai) = nanmedian(z_scores_MSonset_hard);

%         figure(2)
%         scatter(ones(length(z_scores_RT(easyTrials)),1), z_scores_RT(easyTrials))
%         hold on
%         scatter(2*ones(length(z_scores_RT(hardTrials)),1), z_scores_RT(hardTrials))
%         hold on
%         scatter(1, median(z_scores_RT(easyTrials)), 150)
%         hold on
%         scatter(2, median(z_scores_RT(hardTrials)), 150)
%         hold off

        [fnorm_easy, ~] = ksdensity(z_scores_MSonset_easy, x_values, 'function','cdf');
        [fnorm_hard, xinorm] = ksdensity(z_scores_MSonset_hard, x_values, 'function','cdf');

        %total_density_sum = sum(fnorm_easy) + sum(fnorm_hard); %% this makes each probability relative to total events
        f1_normalized = fnorm_easy * (numel(z_scores_MSonset_easy)/possibleEasy);
        f2_normalized = fnorm_hard * (numel(z_scores_MSonset_hard)/possibleHard);

        probDensity{ai}(si,1:length(x_values),1) = f1_normalized; 
        probDensity{ai}(si,1:length(x_values),2) = f2_normalized;

        totalQualifiedTrials(si, ai) = length(postOnsetMS_cond(easyTrials))+length(postOnsetMS_cond(hardTrials));

    end

    if plotOn ~= 0
        sgtitle(subj)
        f1.Position = [1282 335 1184 742];
    end
end


%% Flatten grant array for easy/hard

normalizeData = 0; % this will make the probability distributions to mean "how probable, out of all trials of this trial type, would MS occur at or before X

probDelay = nan(length(subjects), 3, length(analyses)); % 3 time measures saved

% for shaded patch
x_patch = [0, 500, 500, 0]; % x-coordinates
y_patch = [0, 0, 1, 1]; % y-coordinates

%%
for si=1:length(subjects)
    subjectID = si
    subj = subjects{si};
    figure
    figName = sprintf('%s/%s/ProcessedData/Summary/figures/%s', datadir, subj, sprintf('%sprobDensityFunction', subj));
    for ai=1:2 %length(analyses)
    
        analysis_type = analyses{ai};
    
        if strcmp(analysis_type, 'direction')
            fieldNames = {'cardinal', 'oblique'};     
            color = {[17, 119, 51],[51, 34, 136]}; 
        elseif strcmp(analysis_type, 'tilt')
            fieldNames = {'largeoffset','smalloffset'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'location')
            fieldNames = {'horizontalLoc','verticalLoc'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'dirtilt')
            fieldNames = {'easycardinal','hardoblique'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'outcome')
            fieldNames = {'correct','incorrect'};
            color = {[0, 255, 0],[255, 0, 0]}; 
        end
    
        easy_flat = grand_z_scores_MSonset{ai}(subjectID,:,1);
        hard_flat = grand_z_scores_MSonset{ai}(subjectID,:,2);
        
        easy_flat = easy_flat(:); easy_flat = easy_flat(~isnan(easy_flat));
        hard_flat = hard_flat(:); hard_flat = hard_flat(~isnan(hard_flat));
        
        possibleEasy = sum(grand_possible{ai}(subjectID,:,1));
        possibleHard = sum(grand_possible{ai}(subjectID,:,2));
        
        % randomly draw the # of values in each array with replacement
        % fit a kdensity function 1000 times
        x_values = linspace(0, 1500, 500); %linspace(-5, 5, 100);
        numBootstraps = 1000; %500;
        btSamples_easy = nan(numBootstraps, numel(x_values));
        btSamples_hard = nan(numBootstraps, numel(x_values));
        
        for bi=1:numBootstraps
        
        if normalizeData
            ratio2norm_easy = (numel(easy_flat))/possibleEasy;
            ratio2norm_hard = (numel(hard_flat))/possibleHard;
        else
            ratio2norm_easy = 1;
            ratio2norm_hard = 1;
        end
    
            if bi==1
                % mean values (only compute once)
                [fnorm_easy, ~] = ksdensity(easy_flat, x_values, 'function','cdf');
                [fnorm_hard, xinorm] = ksdensity(hard_flat, x_values, 'function','cdf');
                
                 f1_normalized = fnorm_easy * ratio2norm_easy; % this would tell how probably it is across exp
                 f2_normalized = fnorm_hard * ratio2norm_hard;
            end
        
            % bootstrap
            easy_flat_bt = datasample(easy_flat, numel(easy_flat), 'Replace', true);
            hard_flat_bt = datasample(hard_flat, numel(hard_flat), 'Replace', true);
            
            [fnorm_easy, ~] = ksdensity(easy_flat_bt, x_values, 'function','cdf');
            [fnorm_hard, xinorm] = ksdensity(hard_flat_bt, x_values, 'function','cdf');
            
             %f1_normalized 
             btSamples_easy(bi,1:numel(x_values)) = fnorm_easy * ratio2norm_easy; 
             %f2_normalized
             btSamples_hard(bi,1:numel(x_values)) = fnorm_hard * ratio2norm_hard;
        
        end
        
        % Compute the lower and upper percentiles for the 68% confidence interval
        lbd_easy = prctile(btSamples_easy, 2.5,1); %16, 1);
        ubd_easy = prctile(btSamples_easy, 97.5,1); %84, 1);
        lbd_hard = prctile(btSamples_hard, 2.5,1); %16, 1);
        ubd_hard = prctile(btSamples_hard, 97.5,1); %84, 1);
        
         subplot(1,2,ai)
         % Plot the square with a shaded fill
         fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
         hold on
         plot(xinorm, f1_normalized, 'Color', color{1}/255, 'LineWidth',1.5); % Plot the KDE
         hold on
         %errorbar(xinorm, f1_normalized, f1_normalized - lbd_easy, ubd_easy - f1_normalized, 'b', 'LineWidth', 1.5, 'LineStyle', 'none');
         %hold on
         shadedErrorBar(xinorm, f1_normalized, [ubd_easy-f1_normalized; f1_normalized-lbd_easy], 'lineprops', {'-','Color',color{1}/255},'transparent',1,'patchSaturation',0.1);
         hold on
         [half_val_1, half_index_1] = min(abs(f1_normalized - .5));
         time1 = xinorm(half_index_1);
         plot([time1, time1], [0, f1_normalized(half_index_1)], ':', 'Color',color{1}/255, 'LineWidth',4);
         hold on
         %plot([0, time1], [f1_normalized(half_index_1), f1_normalized(half_index_1)], ':', 'Color',color{1}/255, 'LineWidth',2);
         %hold on
         plot(xinorm, f2_normalized, 'Color', color{2}/255, 'LineWidth',1.5); % Plot the KDE
         hold on
         %errorbar(xinorm, f2_normalized, f2_normalized - lbd_hard, ubd_hard - f2_normalized, 'b', 'LineWidth', 1.5, 'LineStyle', 'none');
         %hold on
         shadedErrorBar(xinorm, f2_normalized, [ubd_hard-f2_normalized; f2_normalized-lbd_hard], 'lineprops', {'-','Color',color{2}/255},'transparent',1,'patchSaturation',0.2);
         hold on
         [half_val_2, half_index_2] = min(abs(f2_normalized - .5));
         time2 = xinorm(half_index_2);
         plot([time2, time2], [0, f2_normalized(half_index_2)], ':', 'Color',color{2}/255, 'LineWidth',4);
         hold on
         plot([0, time2], [f2_normalized(half_index_2), f2_normalized(half_index_2)], ':', 'Color','k', 'LineWidth',4);
         %title(analyses{ai})
         %xlim([-3 3])
         probDelay(si,1,ai) = time1;
         probDelay(si,2,ai) = time2;
         probDelay(si,3,ai) = time2-time1;
         %set(gca, 'FontName', 'Arial', 'FontSize', 12);
         
         dpi = get(0, 'ScreenPixelsPerInch');
         set(gca, 'FontName', 'Arial', 'FontSize', (365/dpi)*4.5);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width

         ylabel('cumulative probability density')
         xlabel('time (ms)')
         axis square
         box off
    end

    f = gcf;
    set(f, 'Position', [94 972 893 365])
%     set(f, 'PaperOrientation', 'landscape');

    % Get the current figure position in pixels
%     fig_position_pixels = f.Position;
%     screen_ppi = get(0, 'ScreenPixelsPerInch');
%     fig_position_inches = fig_position_pixels(3:4) / screen_ppi;
%     set(gca, 'FontName', 'Arial', 'FontSize', 12);
%     set(f, 'PaperPosition', [0 0 fig_position_inches/2]);
%         
%     saveas(gca, strcat(figName, '.fig')); % Save as FIG format
%     saveas(gca, strcat(figName, '.pdf')); % Save as PDF format

    
    width_inch = 893 / dpi;
    height_inch = 365 / dpi;

    %pathtoSave = '/Users/rje257/Desktop/oculomotor_manuscript_figures/';
    %set(gcf, 'PaperOrientation', 'landscape');
    %set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    %print(gcf, '-dpdf', fullfile(pathtoSave,'fig6.pdf'));  % for PDF
end

%%

ai= 1;
[nu, val, ci, stats] = ttest(probDelay(:,2,ai), probDelay(:,1,ai))
d = computeCohen_d(probDelay(:,2,ai), probDelay(:,1,ai), 'paired')

ai= 2;
[nu, val, ci, stats] = ttest(probDelay(:,2,ai), probDelay(:,1,ai))
d = computeCohen_d(probDelay(:,2,ai), probDelay(:,1,ai), 'paired')

%% RT vs Probability Density

figure

figName = sprintf('%s/onset_by_rt', datadir);


for ai=1:2 %length(analyses)

    analysis_type = analyses{ai};

    if strcmp(analysis_type, 'direction')
        fieldNames = {'cardinal', 'oblique'};
        color = {[17, 119, 51],[51, 34, 136]}; 
    elseif strcmp(analysis_type, 'tilt')
        fieldNames = {'largeoffset','smalloffset'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'location')
        fieldNames = {'horizontalLoc','verticalLoc'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'dirtilt')
        fieldNames = {'easycardinal','hardoblique'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'outcome')
        fieldNames = {'correct','incorrect'};
        color = {[0, 255, 0],[255, 0, 0]}; 
    end


    %subplot(2,3,ai)
    subplot(1,2,ai)
    
    %[x_fit, y_fit] = fitLine(rtConds(:,1,ai), probDelay(:,1,ai)');
    %plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', color{1}/255);

    %[x_fit2, y_fit2] = fitLine(rtConds(:,2,ai), probDelay(:,2,ai)');
    [x_fit2, y_fit2] = fitLine([rtConds(:,2,ai); rtConds(:,1,ai)], [probDelay(:,2,ai)', probDelay(:,1,ai)']);
    [correlation_coefficient1, pval1] = corr([rtConds(:,2,ai); rtConds(:,1,ai)], [probDelay(:,2,ai); probDelay(:,1,ai)], 'type', 'Spearman', 'rows', 'complete')
    plot(x_fit2, y_fit2, '-', 'LineWidth', 2, 'Color', 'k'); %color{2}/255);

    hold on
    %scatter(rtConds(:,3,ai), probDelay(:,3,ai), 'filled', 'MarkerFaceColor', 'k')
    scatter(rtConds(:,1,ai), probDelay(:,1,ai), 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
    hold on
    scatter(rtConds(:,2,ai), probDelay(:,2,ai), 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
    hold on
    scatter(mean(rtConds(:,1,ai)), mean(probDelay(:,1,ai)), 500, '+', 'MarkerEdgeColor', color{1}/255, 'LineWidth', 5)
    hold on
    scatter(mean(rtConds(:,2,ai)), mean(probDelay(:,2,ai)), 500, '+', 'MarkerEdgeColor', color{2}/255, 'LineWidth', 5)

    hold on

    % plot the identity line
    x_values = linspace(0, 10, 100);
    y_values = x_values;
    plot(x_values, y_values, 'k--');
    axis equal
    set(gca, 'LineWidth', 2);
    %grid on
    xlim([500 1000])
    ylim([300 800])
    ax1 = gca;
    ax1.YTick = 300:100:800;
    ax1.XTick = 500:100:1000;
    %set(gca, 'FontName', 'Arial', 'FontSize', 20);
    box off
    xlabel('response time (ms)')
    ylabel('microsaccade rebound (ms)')
    set(gca, 'LineWidth', 2);

    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*1.8);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width

%     subplot(2,3,ai+3)
% 
%     [x_fit, y_fit] = fitLine(rtConds(:,3,ai), probDelay(:,3,ai)');
%     %correlation_coefficient1 = corr(rtConds(:,3,ai), probDelay(:,3,ai), 'type', 'Spearman', 'rows', 'complete')
%     plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'b');
%     hold on
%     scatter(rtConds(:,3,ai), probDelay(:,3,ai), 85, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w','LineWidth', 1, 'MarkerFaceAlpha', 0.5)
%     xlabel('response time (hard - easy) (ms)')
%     ylabel('estimated time of first msacc (hard - easy)')
%     set(gca, 'FontName', 'Arial', 'FontSize', 12);
%     axis equal
%     %grid on
%     ylim([-40 160])
%     xlim([-50 150])
%     ax2 = gca;
%     box off
%     ax2.YTick = -40:50:140;
%     set(gca, 'FontName', 'Arial', 'FontSize', 12);
end

% f = gcf;
% set(f, 'Position', [848 426 1672 911])
% 
% set(f, 'PaperOrientation', 'landscape');

% Get the current figure position in pixels
%fig_position_pixels = f.Position;
%screen_ppi = get(0, 'ScreenPixelsPerInch');
%fig_position_inches = fig_position_pixels(3:4) / screen_ppi;
%set(f, 'PaperPosition', [0 0 fig_position_inches/2]);
   

f1 = gcf;
f1.Position = [334 778 889 513];

pathtoSave = '/Users/rje257/Desktop/oculomotor_manuscript_figures/';
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gcf, 'PaperPositionMode', 'auto'); 

%set(gca, 'LineWidth', 1);
print(gcf, '-dpdf', fullfile(pathtoSave,'fig7.pdf'));  % for PDF

% saveas(gca, strcat(figName, '.fig')); % Save as FIG format
% saveas(gca, strcat(figName, '.pdf')); % Save as PDF format

%% Pupil vs Probability Delay (MS)

figure

event = 2;
paramType = 'peak'; %'peak'; % peak or amp

figName = sprintf('%s/onset_by_%speak%s', datadir, paramType, num2str(event));

for ai=1:2 %length(analyses)

    analysis_type = analyses{ai};

    if strcmp(analysis_type, 'direction')
        fieldNames = {'cardinal', 'oblique'};
        color = {[17, 119, 51],[51, 34, 136]}; 
    elseif strcmp(analysis_type, 'tilt')
        fieldNames = {'largeoffset','smalloffset'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'location')
        fieldNames = {'horizontalLoc','verticalLoc'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'dirtilt')
        fieldNames = {'easycardinal','hardoblique'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'outcome')
        fieldNames = {'correct','incorrect'};
        color = {[0, 255, 0],[255, 0, 0]}; 
    end


    %subplot(2,3,ai)
    subplot(1,2,ai)

    if strcmp(paramType, 'peak')
        peak1 = (latvals{ai}(:,event,1)+eventtimes{ai}(:,event,1)+930);
        peak2 = (latvals{ai}(:,event,2)+eventtimes{ai}(:,event,2)+930);
    elseif strcmp(paramType, 'amp')
        peak1 = (amp{ai}(:,event,1));
        peak2 = (amp{ai}(:,event,2));
    elseif strcmp(paramType, 'boxamp')
        peak1 = (boxamp{ai}(:,event,1));
        peak2 = (boxamp{ai}(:,event,2));
    end

    [x_fit2, y_fit2] = fitLine([peak2; peak1], [probDelay(:,2,ai)', probDelay(:,1,ai)']);

    plot(x_fit2, y_fit2, '-', 'LineWidth', 2, 'Color', 'k'); %color{2}/255);

    hold on
    scatter(peak1, probDelay(:,1,ai), 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
    hold on
    scatter(peak2, probDelay(:,2,ai), 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
    hold on
    scatter(mean(peak1), mean(probDelay(:,1,ai)), 500, '+', 'MarkerEdgeColor', color{1}/255, 'LineWidth', 5)
    hold on
    scatter(mean(peak2), mean(probDelay(:,2,ai)), 500, '+', 'MarkerEdgeColor', color{2}/255, 'LineWidth', 5)
    [correlation_coefficient1, pval1] = corr([peak1; peak2], [probDelay(:,1,ai); probDelay(:,2,ai)], 'type', 'Spearman', 'rows', 'complete')
    hold on

    % plot the identity line
    x_values = linspace(0, 10, 100);
    y_values = x_values;
    plot(x_values, y_values, 'k--');
    %grid on
    set(gca, 'FontName', 'Arial', 'FontSize', 20);
    ax1 = gca;
    if strcmp(paramType, 'peak')
        xlabel('peak dilation (ms)', 'FontSize', 20)
        if event ==1
            xlim([150 500])
            ax1.XTick = 150:200:500;
        elseif event==2
            xlim([1050 1850])
            ylim([150 950])
            ax1.XTick = 1050:200:1850;
        end
        ylim([300 800])
        ax1.YTick = 200:200:800;
        
    elseif strcmp(paramType, 'amp')
        xlabel('amp (psc)')
        xlim([0 150])
        ax1.XTick = 0:40:150;
         ylim([150 950])
        ax1.YTick = 150:200:950;
    elseif strcmp(paramType, 'boxamp')
        xlabel('boxamp (psc)')
        xlim([-20 100])
        ax1.XTick = -20:40:100;
         ylim([200 900])
        ax1.YTick = 200:200:900;
    end
    ylabel('microsaccade rebound (ms)')
    box off
    set(gca, 'LineWidth', 2);

    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*1.8);
    ax = gca;  % Get the current axis
    ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
    ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width

%     subplot(2,3,ai+3)
% 
%     diff_peak = peak2 - peak1;
%     [x_fit, y_fit] = fitLine(diff_peak, probDelay(:,3,ai)');
%     %correlation_coefficient1 = corr(diff_peak, probDelay(:,3,ai), 'type', 'Spearman', 'rows', 'complete')
%     plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'b');
%     hold on
%     scatter(diff_peak, probDelay(:,3,ai), 85, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w','LineWidth', 1, 'MarkerFaceAlpha', 0.5)
%     ax2 = gca;
%     if strcmp(paramType, 'peak')
%         if event ==1
%             xlim([-175 100])
%         elseif event ==2
%             xlim([-25 200])
%         elseif event == 3
%             xlim([-200 150])
%         end
%         xlabel('peak time (hard - easy) (ms)')
%         ylim([-65 160])
%         ax2.YTick = -65:50:160;
%     elseif strcmp(paramType, 'amp')
%         xlabel('amp (hard - easy) (psc)')
%         xlim([-10 25])
%     elseif strcmp(paramType, 'boxamp')
%         xlabel('boxamp (hard - easy) (psc)')
%         xlim([-5 35])
%     end
%     ylabel('estimated time of first msacc (hard - easy)')
    %grid on
    axis square

end

% f = gcf;
% set(f, 'Position', [848 426 1672 911])
% 
% 
% set(f, 'PaperOrientation', 'landscape');
% 
% % Get the current figure position in pixels
% fig_position_pixels = f.Position;
% screen_ppi = get(0, 'ScreenPixelsPerInch');
% fig_position_inches = fig_position_pixels(3:4) / screen_ppi;
% set(f, 'PaperPosition', [0 0 fig_position_inches/2]);

% set(gca, 'FontName', 'Arial', 'FontSize', 20);
% set(gca, 'LineWidth', 2);
% 
% f1 = gcf;
% f1.Position = [128 981 1020 309];
% 
% saveas(gca, strcat(figName, '.fig')); % Save as FIG format
% saveas(gca, strcat(figName, '.pdf')); % Save as PDF format


f1 = gcf;
f1.Position = [334 778 889 513];

pathtoSave = '/Users/rje257/Desktop/oculomotor_manuscript_figures/';
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gcf, 'PaperPositionMode', 'auto'); 

%set(gca, 'LineWidth', 1);
print(gcf, '-dpdf', fullfile(pathtoSave,'fig11.pdf'));  % for PDF


%%
figure(1)

n = length(subjects);
ones_vector = zeros(n, 1);
jittered_vector = ones_vector + (rand(n, 1) - 0.5) * 0.55;

figName1 = sprintf('%s/ampPupilperCond', datadir);
figName2 = sprintf('%s/peakPupilperCond', datadir);
figName3 = sprintf('%s/ampDiffperCond', datadir);
figName4 = sprintf('%s/peakDiffperCond', datadir);

for ai=1:2 %2 %length(analyses)

    analysis_type = analyses{ai};

    if strcmp(analysis_type, 'direction')
        fieldNames = {'cardinal', 'oblique'};
        color = {[17, 119, 51],[51, 34, 136]}; 
        alphlev1 = 0.5;
        lincol1 = 'k';
        alphlev2 = 0.5;
        lincol2 = 'k';
        condComparison = '(obl - card)';
    elseif strcmp(analysis_type, 'tilt')
        fieldNames = {'largeoffset','smalloffset'};
        color = {[0, 0, 0],[175, 175, 175]}; 
        alphlev1 = 0.5;
        lincol1 = 'k';
        alphlev2 = 0.5;
        lincol2 = 'k';
        condComparison = '(small - large)';
    elseif strcmp(analysis_type, 'location')
        fieldNames = {'horizontalLoc','verticalLoc'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'dirtilt')
        fieldNames = {'easycardinal','hardoblique'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'outcome')
        fieldNames = {'correct','incorrect'};
        color = {[0, 255, 0],[255, 0, 0]}; 
    end

%     figure(1)
%     %subplot(1,length(analyses),ai)
%     subplot(1,2,ai)
%     plot([1 2], [amp{ai}(:,1,1), amp{ai}(:,1,2)], 'k', 'LineWidth', 1)
%     hold on
%     scatter(ones(length(amp{ai}(:,1,1)),1), amp{ai}(:,1,1), 200, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     hold on
%     scatter(2*ones(length(amp{ai}(:,1,2)),1), amp{ai}(:,1,2), 200, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1)
%     hold on
%     plot([4 5], [amp{ai}(:,2,1), amp{ai}(:,2,2)], 'k', 'LineWidth', 1)
%     hold on
%     scatter(4*ones(length(amp{ai}(:,2,1)),1), amp{ai}(:,2,1), 200, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     hold on
%     scatter(5*ones(length(amp{ai}(:,2,2)),1), amp{ai}(:,2,2), 200, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     hold on
%     plot([7 8], [amp{ai}(:,3,1), amp{ai}(:,3,2)], 'k', 'LineWidth', 1)
%     hold on
%     scatter(7*ones(length(amp{ai}(:,3,1)),1), amp{ai}(:,3,1), 200, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     hold on
%     scatter(8*ones(length(amp{ai}(:,3,2)),1), amp{ai}(:,3,2), 200, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     hold on
%     plot([10 11], [boxamp{ai}(:,1,1), boxamp{ai}(:,1,2)], 'k', 'LineWidth', 1)
%     hold on
%     scatter(10*ones(length(boxamp{ai}(:,1,1)),1), boxamp{ai}(:,1,1), 200, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     hold on
%     scatter(11*ones(length(boxamp{ai}(:,1,2)),1), boxamp{ai}(:,1,2), 200, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 1) 
%     xlim([0,12])
%     ax1 = gca;
%     ax1.XTick = linspace(1.5, 10.5, 4);
%     ylabel('amplitude (psc)')
%     ax1.XTickLabel = {"trial start"; "task-evoked"; "response-evoked"; "stimulus boxcar"};
%     f1 = gcf;
%     set(gca, 'FontName', 'Arial', 'FontSize', 20);
%     set(gca, 'LineWidth', 2);
%     box off

    figure(2)
    %subplot(1,length(analyses),ai)
    subplot(1,2,ai)
    [nu, val, ci, stats] = ttest((latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930)-(latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930))
    %d = computeCohen_d((latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930), (latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930), 'paired')
    hold on
    scatter(jittered_vector+ones(length(latvals{ai}(:,1,1)),1), latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930, 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', lincol1, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev1)
    hold on
    scatter(jittered_vector+(3*ones(length(latvals{ai}(:,1,2)),1)), latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930, 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', lincol2, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev2)
    hold on
    plot([jittered_vector+1 jittered_vector+3]', [latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930, latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930]', 'k', 'LineWidth', 1)
    hold on
    [nu, val, ci, stats] = ttest((latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930)-(latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930))
    %d = computeCohen_d((latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930), (latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930), 'paired')
    hold on
    scatter(jittered_vector+(5*ones(length(latvals{ai}(:,2,1)),1)), latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930, 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', lincol1, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev1)
    hold on
    scatter(jittered_vector+(7*ones(length(latvals{ai}(:,2,2)),1)), latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930, 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', lincol2, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev2)
    hold on
    plot([jittered_vector+5 jittered_vector+7]', [latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930, latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930]', 'k', 'LineWidth', 1)
    hold on
    [nu, val, ci, stats] = ttest((latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930)-(latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930))
     %d = computeCohen_d((latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930), (latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930), 'paired')
    scatter(jittered_vector+(9*ones(length(latvals{ai}(:,3,1)),1)), latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930, 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', lincol1, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev1)
    hold on
    scatter(jittered_vector+(11*ones(length(latvals{ai}(:,3,2)),1)), latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930, 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', lincol2, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev2)
    hold on
    plot([jittered_vector+9 jittered_vector+11]', [latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930, latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930]', 'k', 'LineWidth', 1)
    ax1 = gca;
    ax1.XTick = linspace(2, 10, 3);
    ylabel('peak dilation (ms)')
    ax1.XTickLabel = {"trial start"; "task-evoked"; "response-evoked"};
    xlim([0,12])
    f2 = gcf;
    set(gca, 'XTickLabelRotation', 45);
    %set(gca, 'FontName', 'Arial', 'FontSize', 20);
    set(gca, 'LineWidth', 2);
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*1.8);
    ax = gca;  % Get the current axis
    ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
    ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    box off
    axis square

    figure(3)
    %subplot(1,length(analyses),ai)
    subplot(1,2,ai)
    [nsubs, ~] = size(amp{ai}(:,1,1));
    psc1 = ((amp{ai}(:,1,2).*psc_scaleFactor)-(amp{ai}(:,1,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc1)
    d = computeCohen_d(amp{ai}(:,1,2),amp{ai}(:,1,1), 'paired')
    hold on
    sed1 = nanstd(psc1) / sqrt(nsubs);
    scatter(jittered_vector+(1.5*ones(length(amp{ai}(:,1,1)),1)), psc1, 250, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2) 
    hold on
    errorbar(1.5, nanmean(psc1), sed1, 'k-', 'LineWidth',5, 'CapSize', 15);
    hold on
    %scatter(1.5, nanmean(psc1), 200, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 3);
    hold on
    %boxplot(amp{ai}(:,1,2)-amp{ai}(:,1,1), 'Whisker', 1.5, 'Positions', 1.5, 'Widths', 1);
    hold on
    psc2 = ((amp{ai}(:,2,2).*psc_scaleFactor)-(amp{ai}(:,2,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc2)
    d = computeCohen_d(amp{ai}(:,2,2),amp{ai}(:,2,1), 'paired')
    sed2 = nanstd(psc2) / sqrt(nsubs);
    hold on
    scatter(jittered_vector+(5.5*ones(length(amp{ai}(:,2,1)),1)), psc2, 250, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2) 
    hold on
    errorbar(5.5, nanmean(psc2), sed2, 'k-', 'LineWidth',5, 'CapSize', 15);
    hold on
    %scatter(4.5, nanmean(psc2), 200, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 3);
    hold on
    %boxplot(amp{ai}(:,2,2)-amp{ai}(:,2,1), 'Whisker', 1.5, 'Positions', 4.5, 'Widths', 1);
    hold on   
    psc3 = ((amp{ai}(:,3,2).*psc_scaleFactor)-(amp{ai}(:,3,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc3)
    d = computeCohen_d(amp{ai}(:,3,2),amp{ai}(:,3,1), 'paired')
    sed3 = nanstd(psc3) / sqrt(nsubs);
    scatter(jittered_vector+(9.5*ones(length(amp{ai}(:,3,1)),1)), psc3, 250, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2) 
    hold on
    errorbar(9.5, nanmean(psc3), sed3, 'k-', 'LineWidth',5, 'CapSize', 15);
    hold on
    %scatter(7.5, nanmean(psc3), 200, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 3);
    hold on
    %boxplot(amp{ai}(:,3,2)-amp{ai}(:,3,1), 'Whisker', 1.5, 'Positions', 7.5, 'Widths', 1);
    psc4 = ((boxamp{ai}(:,1,2).*psc_scaleFactor)-(boxamp{ai}(:,1,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc4)
    d = computeCohen_d(boxamp{ai}(:,1,2),boxamp{ai}(:,1,1), 'paired')
    sed4 = nanstd(psc4) / sqrt(nsubs);
    scatter(jittered_vector+(13.5*ones(length(boxamp{ai}(:,1,1)),1)), psc4, 250, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2)
    hold on
    errorbar(13.5, nanmean(psc4), sed4, 'k-', 'LineWidth',5, 'CapSize', 15);
    hold on 
    %scatter(10.5, nanmean(psc4), 200, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 3);
    hold on
    %boxplot(boxamp{ai}(:,1,2)-boxamp{ai}(:,1,1), 'Whisker', 1.5, 'Positions', 10.5, 'Widths', 1);
    hold on
    yline(0, ':k', 'LineWidth',2)
    xlim([0,15])
    ylim([-15,35])
    ax1 = gca;
    ax1.XTick = linspace(1.5, 13.5, 4);
    ylabel(['pupil area % ', condComparison]);
    set(gca, 'XTickLabelRotation', 45);
    ax1.XTickLabel = {"trial start"; "task-evoked"; "response-evoked"; "stimulus boxcar"};
    f3 = gcf;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*1.8);
    ax = gca;  % Get the current axis
    ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
    ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    box off
    %axis square

%     figure(4) % PEAK DIFF
%     subplot(1,length(analyses),ai)
%     [nsubs, ~] = size(amp{ai}(:,1,1));
%     psc1 = (latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930) - (latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930);
%     hold on
%     sed1 = nanstd(psc1) / sqrt(nsubs);
%     scatter(1.5*ones(length(latvals{ai}(:,1,1)),1), psc1, 85, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 1, 'MarkerFaceAlpha', 0.2) 
%     hold on
%     errorbar(1.5, nanmean(psc1), sed1, 'k-', 'LineWidth',3, 'CapSize', 15);
%     hold on
%     psc2 = (latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930) - (latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930);
%     scatter(4.5*ones(length(latvals{ai}(:,2,1)),1), psc2, 85, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 1, 'MarkerFaceAlpha', 0.2) 
%     hold on
%     sed2 = nanstd(psc2) / sqrt(nsubs);
%     errorbar(4.5, nanmean(psc2), sed2, 'k-', 'LineWidth',3, 'CapSize', 15);
%     hold on
%     %boxplot(amp{ai}(:,2,2)-amp{ai}(:,2,1), 'Whisker', 1.5, 'Positions', 4.5, 'Widths', 1);
%     hold on   
%     psc3 = (latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930) - (latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930);
%     sed3 = nanstd(psc3) / sqrt(nsubs);
%     scatter(7.5*ones(length(amp{ai}(:,3,1)),1), psc3, 85, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 1, 'MarkerFaceAlpha', 0.2) 
%     hold on
%     errorbar(7.5, nanmean(psc3), sed3, 'k-', 'LineWidth',3, 'CapSize', 15);
%     yline(0, ':k', 'LineWidth',2)
%     xlim([0,9])
%     ax1 = gca;
%     ax1.XTick = linspace(1.5, 7.5, 3);
%     ylabel('peak difference (ms)')
%     ax1.XTickLabel = {"trial start"; "task-evoked"; "response-evoked"};
%     f4 = gcf;
%     set(gca, 'FontName', 'Arial', 'FontSize', 20);
%     set(gca, 'LineWidth', 2);
%     box off

end

for fi=[f1, f2,f3] %, f4]

    %f1 = gcf;
    %f1.Position = [334 778 889 513];
    if fi==1
        figName = 'x'; %figName1;
        f = figure(1);
        set(f, 'Position', [334 778 889 513]); %[199 465 1458 872])
    elseif fi ==2
        figName = '9'; figName2;
        f = figure(2);
        set(f, 'Position', [334 778 889 513]); %[190 632 1639 648])
    elseif fi ==3
        figName = '10'; %figName3;
        f = figure(3);
        set(f, 'Position', [194 784 956 426]); %[67 547 1971 739])
    elseif fi ==4
        figName = 'x'; %figName4;
        f = figure(4);
    end

%     
%     set(f, 'PaperOrientation', 'landscape');
%     
%     % Get the current figure position in pixels
%     fig_position_pixels = f.Position;
%     screen_ppi = get(0, 'ScreenPixelsPerInch');
%     fig_position_inches = fig_position_pixels(3:4) / screen_ppi;
%     set(f, 'PaperPosition', [0 0 fig_position_inches/2]);
        
%     set(gca, 'FontName', 'Arial', 'FontSize', 20);
%     set(gca, 'LineWidth', 2);

    %f1 = gcf;
    %f1.Position = [334 778 889 513]; %[128 888 844 402];

    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    set(gcf, 'PaperPositionMode', 'auto'); 

    %set(gca, 'LineWidth', 1);
    tmp = strsplit(figName, '/');
    %print(gcf, '-dpdf', fullfile(pathtoSave,sprintf('fig%s.pdf', figName)));  % for PDF tmp{8}))); %

    %saveas(gca, strcat(figName, '.fig')); % Save as FIG format
    %saveas(gca, strcat(figName, '.pdf')); % Save as PDF format
end

pathtoSave = '/Users/rje257/Desktop/oculomotor_manuscript_figures/';
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
%set(gcf, 'PaperPositionMode', 'auto'); 

%set(gca, 'LineWidth', 1);
print(gcf, '-dpdf', fullfile(pathtoSave,'fig9.pdf'));  % for PDF

%% raw data

figure
subplot(1,2,1)
scatter(ones(length(msOnset(:,1,1))), msOnset(:,1,1), 'g')
hold on
scatter(2*ones(length(msOnset(:,2,1))), msOnset(:,2,1), 'r')
hold on
plot([msOnset(:,1,1), msOnset(:,2,1)]')
xlim([0 3])
subplot(1,2,2)
scatter(ones(length(msOnset(:,1,2))), msOnset(:,1,2), 'g')
hold on
scatter(2*ones(length(msOnset(:,2,2))), msOnset(:,2,2), 'r')
hold on
plot([msOnset(:,1,2), msOnset(:,2,2)]')
xlim([0 3])

%% Difference plot

figure
title('Difference')
scatter(ones(length(msOnset(:,1,2))), msOnset(:,2,1) - msOnset(:,1,1), 'g')
hold on
scatter(2*ones(length(msOnset(:,2,2))), msOnset(:,2,2) - msOnset(:,1,2), 'r')
xlim([0 3])

%% z-scored data (multiply by *1 to make z-score - "early values"
% should not do this for duration

switchSign = 1; % fixed above now

figure
subplot(1,2,1)
scatter(ones(length(zscored_msOnset(:,1,1))), switchSign*zscored_msOnset(:,1,1), 'g')
hold on
scatter(2*ones(length(zscored_msOnset(:,2,1))), switchSign*zscored_msOnset(:,2,1), 'r')
hold on
plot([switchSign*zscored_msOnset(:,1,1), switchSign*zscored_msOnset(:,2,1)]', 'k')
xlim([0 3])
subplot(1,2,2)
scatter(ones(length(zscored_msOnset(:,1,2))), switchSign*zscored_msOnset(:,1,2), 'g')
hold on
scatter(2*ones(length(zscored_msOnset(:,2,2))), switchSign*zscored_msOnset(:,2,2), 'r')
hold on
plot([switchSign*zscored_msOnset(:,1,2), switchSign*zscored_msOnset(:,2,2)]', 'k')
xlim([0 3])

%%

figure

for ai=1:length(analyses)


    analysis_type = analyses{ai};

    % Correalte MS latency with Pupil Latency
    if strcmp(analysis_type, 'direction')
        fieldNames = {'cardinal', 'oblique'};
        color = {[17, 119, 51],[51, 34, 136]}; 
    elseif strcmp(analysis_type, 'tilt')
        fieldNames = {'largeoffset','smalloffset'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'location')
        fieldNames = {'horizontalLoc','verticalLoc'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'dirtilt')
        fieldNames = {'easycardinal','hardoblique'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'outcome')
        fieldNames = {'correct','incorrect'};
        color = {[0, 255, 0],[255, 0, 0]}; 
    end


    currDensity_easy = probDensity{ai}(:,:,1);
    currDensity_hard = probDensity{ai}(:,:,2);

%     if switchSign == -1
%         currDensity_easy = fliplr(currDensity_easy);
%         currDensity_hard = fliplr(currDensity_hard);
%     end
    
    
    grandTotalTrials = sum(totalQualifiedTrials(ai,1));
    
    subjectWeights = totalQualifiedTrials(:,1)./grandTotalTrials;
    %subjectWeights = [1; 0; 0; 0; 0; 0 ; 0; 0];

    % adjust probability density to reflect weighted probability (so I am not
    % placing too much weight on subjects with very few trials)
    currDensity_easy = currDensity_easy.*subjectWeights;
    currDensity_hard = currDensity_hard.*subjectWeights;

    subplot(1,3,ai)
    
    % weighted mean
    easyDensitySeries = sum(currDensity_easy,1); %/(sum(mean(currDensity_easy,1))+sum(mean(currDensity_hard,1)));
    hardDensitySeries = sum(currDensity_hard,1); %/(sum(mean(currDensity_hard,1))+sum(mean(currDensity_easy,1)));

    % Calculate the standard error of the weighted mean (SEWM)
    numerator = sum(subjectWeights .* (currDensity_easy - easyDensitySeries).^2);
    denominator = (sum(subjectWeights))^2;
    sewm_easy = sqrt(numerator / denominator);
    numerator = sum(subjectWeights .* (currDensity_hard - hardDensitySeries).^2);
    sewm_hard = sqrt(numerator / denominator);

    plot(xinorm, easyDensitySeries, 'Color', color{1}/255, 'LineWidth',2); % Plot the KDE
    hold on
    %[~, maxIdx] = max(sum(currDensity_easy,1));
    %xline(median(sum(currDensity_easy,1)), '--', 'Color', color{1}/255, 'LineWidth',2)
    %hold on
    %xline(x_values(maxIdx), 'Color', color{1}/255, 'LineWidth',2)
    %hold on
    plot(xinorm, hardDensitySeries, 'Color', color{2}/255, 'LineWidth',2); % Plot the KDE
    hold on
    %[~, maxIdx] = max(sum(currDensity_hard,1));
    %xline(x_values(maxIdx), 'Color', color{2}/255, 'LineWidth',2)
    %hold on
    %xline(median(sum(currDensity_hard,1)), '--', 'Color', color{2}/255, 'LineWidth',2)
    %hold on
    %errorbar(xinorm, easyDensitySeries, sewm_easy, 'r');
    xlabel('Time relative to stim offset (std normalized)');
    ylabel('probability density');
    title('Kernel Density Estimate');
    title(analyses{ai})
    xlim([-500 500])
end

%%

function [postsupPeakTime, suppressionTime] = calculateSuppPeak(rate, selection, window_size)
    meanMS = nanmean(rate(selection,:),1);

    window_size = 50; % override with 50 for testing

    smoothed_data = smoothdata(meanMS, 'movmean', window_size);
        
    decreasing_index = find(diff(smoothed_data(1800:end)) < 0, 1);
    postsupPeakTime = decreasing_index+1800;

    %[~,suppressionTime] = min(meanMS(1300:1800)); 
    %[~,suppressionTime] = min(smoothed_data(1300:1800)); 
    
    increasing_index = find(diff(smoothed_data(1300:end)) > 0, 1, 'first');
    suppressionTime = increasing_index+1300;

end