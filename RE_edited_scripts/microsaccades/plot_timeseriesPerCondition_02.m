clc; clear all; %close all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
nSampleCutOff = 4000; % or 2500; %
analysis_type = 'direction'; %'direction'; % 'tilt'; %'location'

 % trialwise will take mean over all trials across sessions and bootstrap
 % with 68% CI
 % if 0, will do SEM across sessions and do permutation test
trialwise = 1;
percentileRange = 68; % for plotting, do 68% CIs

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
end

direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
time=[0,100];

x = (1 : nSampleCutOff) - 1300;
samplingRateData = 1000; % 1000 hz after running preprocess_MA_01.m

for ii = 1:1 %8
    figure

    summaryMSPath = fullfile(datadir, subjects{ii}, 'ProcessedData', 'Summary', 'microsaccades');
    summaryFigPath = fullfile(datadir, subjects{ii}, 'ProcessedData', 'Summary', 'figures');

    % for saving figures / variables
    fileName = sprintf('%s_%s_%s_%s_%i_%i', subjects{ii}, analysis_type, fieldNames{1}, fieldNames{2}, trialwise, percentileRange);
    figName = fullfile(summaryFigPath, fileName);
    summaryName = fullfile(summaryMSPath, fileName);

    if isfile(strcat(summaryName, '.mat'))
        load(strcat(summaryName, '.mat'))
        dataReady = 1;
    else
    
        dataReady = 0;
        mkdir(summaryMSPath)
        mkdir(summaryFigPath)
        
        sub_d_car = [];
        sub_d_obl = [];
        RT_marker_car = [];
        RT_marker_obl = [];
        RTcard_all = []; % all RTs including max
        RTobl_all = []; % all RTs including max
        
        rate = struct();
        rateSummary = struct();
        rtSummary = struct();
        rtTrialwise = struct();
    
        for jj = 1 : 2 %2
            if jj== 1 & ii > 6 
                continue
            end
            for kk = 1 : 8
    
                processeddata_folder = fullfile(datadir, subjects{ii}, ...
                    'ProcessedData', protocols{jj});
                MSpath = fullfile(processeddata_folder, 'eyedata','microsaccades');
                MATpath = fullfile(processeddata_folder, 'eyedata','MATs');
    
                dirname = dir(fullfile(MATpath, sprintf('%s*_tab.mat', direction{kk})));
    
                % if I want to include S04 extra data, I should match the tab
                % file
                for di=1:length(dirname) % for multiple blocks % 
    
                    directionFilename = dirname(di).name; % for multiple, it automatically sorts alphanumerically
    
                    if di == 1
                        ms_path= fullfile(MSpath,sprintf('%s_microsaccadeMatrix.mat', direction{kk}));
                    elseif di == 2
                        ms_path= fullfile(MSpath,sprintf('%s2_microsaccadeMatrix.mat', direction{kk}));
                    end

                    try % in one case, eye tracker broke for first session, so only session 2 exists (S03)
                        load(ms_path);
                    catch
                        ms_path= fullfile(MSpath,sprintf('%s2_microsaccadeMatrix.mat', direction{kk}));
                        load(ms_path);
                    end
        
                    tab_path = fullfile(MATpath, directionFilename);
                    load(tab_path)
                    
                    % Make this to filter card-oblique, or by location, or by
                    % tilt (or a combination of these)
    
                    for pp=1:length(fieldNames)
                        fieldName = fieldNames{pp};
    
                        if ~isfield(rateSummary, fieldName)
                            rateSummary.(fieldName) = []; rtSummary.(fieldName) = []; rtTrialwise.(fieldName) = []; ...
                                rate.(fieldName) = []; rtUnfiltered.(fieldName) = []; trialData.(fieldName) = [];
                        end
    
                        if strcmp(fieldName, 'cardinal')
                            locations = 0:45:315; directions = 0:90:270; tilts = [0.5, 1, 2, 4, 8];
                        elseif strcmp(fieldName, 'oblique')
                            locations = 0:45:315; directions = 45:90:315; tilts = [0.5, 1, 2, 4, 8];
                        elseif strcmp(fieldName, 'largeoffset')
                            locations = 0:45:315; directions = 0:45:315; tilts = [8];
                        elseif strcmp(fieldName, 'smalloffset')
                            locations = 0:45:315; directions = 0:45:315; tilts = [0.5, 1, 2, 4];
                        elseif strcmp(fieldName, 'horizontalLoc')
                            locations = [0 180]; directions = 0:45:315; tilts = [8];
                        elseif strcmp(fieldName, 'verticalLoc')
                            locations = [90 270]; directions = 0:45:315; tilts = [0.5, 1, 2, 4];
                        elseif strcmp(fieldName, 'easycardinal')
                            locations = 0:45:315; directions = 0:90:270; tilts = [8];
                        elseif strcmp(fieldName, 'hardoblique')    
                            locations = 0:45:315; directions = 45:90:315; tilts = [0.5, 1, 2, 4];
                        end
    
                        [s_mean,RT_marker, RT_trialwise, RT_unfiltered, binaryRate, tabTrials] = computeMSRateRT(MS_TEMP,tab,samplingRateData, locations, directions, tilts, nSampleCutOff);
    
                        if ~all(isnan(s_mean(1000:end))) % checks if the condition is probed in this loop
                            rate.(fieldName) = [rate.(fieldName);binaryRate(:,1000:end)]; % all trials concatenated
                            rateSummary.(fieldName) = [rateSummary.(fieldName);s_mean(1000:end)]; % summary per session
                            rtSummary.(fieldName) = [rtSummary.(fieldName), RT_marker]; 
                            rtTrialwise.(fieldName) = [rtTrialwise.(fieldName); RT_trialwise]; % *1000 to convert ms to s 
                            rtUnfiltered.(fieldName) = [rtUnfiltered.(fieldName); RT_unfiltered]; % all RTs (no filter)
                            trialData.(fieldName) = [trialData.(fieldName); tabTrials];
                        end
                    end
                      
                end
            end
        end
    end

    %% subjectwise temporal plot

    for pp=1:length(fieldNames)
        fieldName = fieldNames{pp};

        % Plotting the confidence intervals (optional)
        if trialwise
            bootFile = strcat(summaryName, sprintf('%sbootstraps.mat', fieldName));
            if isfile(bootFile)
                load(bootFile)
            else
                disp('Bootstrapping Data .. ')
                bootstrapStatistics = bootstrapData(rate, fieldName);
                save(strcat(summaryName, sprintf('%sbootstraps.mat', fieldName)), 'bootstrapStatistics');
            end

            [lowerBound, upperBound] = findCI(bootstrapStatistics, percentileRange);
            upperError = upperBound-nanmean(rate.(fieldName),1)*60;
            lowerError = nanmean(rate.(fieldName),1)*60-lowerBound;
            
            a = shadedErrorBar(x(1000:end), nanmean(rate.(fieldName),1)*60, [upperError;lowerError], 'lineprops', {'-','Color',color{pp}/255},'transparent',1,'patchSaturation',0.2);
            hold on
            b = plot(x(1000:end), nanmean(rate.(fieldName),1)*60,'-','LineWidth',1, 'Color', color{pp}/255);
            hold on
            xline(nanmean(rtTrialwise.(fieldName)),':', 'Mean RT', 'Color', color{pp}/255, 'LineWidth',2) % all trials concatenated
        
        % Plotting SEM across runs
        else
            semVal = std(rateSummary.(fieldName),0,1)/sqrt(size(rateSummary.(fieldName),1));
            a = shadedErrorBar(x(1000:end), mean(rateSummary.(fieldName),1),semVal, 'lineprops', {'-','Color',color{pp}/255},'transparent',1,'patchSaturation',0.2);
            hold on
            b = plot(x(1000:end), mean(rateSummary.(fieldName),1),'-','LineWidth',1, 'Color', color{pp}/255);
            hold on
            xline(mean(rtSummary.(fieldName)),':', 'Mean RT', 'Color', color{pp}/255, 'LineWidth',2) % session means
        end

        % plot trial STIM on and STIM off
        hold on
        xline(0,'-','Stimuli on')
        hold on 
        xline(500,'-','Stimuli off')
        hold on

%         b.Color = color{pp}/255;
        ylim([0,6])
        xlim([-300,2000])
        ylabel('microsaccade rate (hz)')
        xlabel('time (ms)')
        title(sprintf('%s %s %s vs. %s', subjects{ii}, analysis_type, fieldNames{1}, fieldNames{2}))
        f = gcf;
        f.Position = [5 996 1431 311];

        yLimits = ylim; % get current ylim

    end


    % plot significant clusters (this method only when doing trialwise)

    if length(fieldNames) == 2 %&& ~trialwise % do cluster permutation
        sig_cluster = checksignificance_perm(rateSummary);
    %elseif length(fieldNames) == 2 && trialwise % check for overlap but correct for MC
    end

    for si=1:length(sig_cluster)
        segLength = length(sig_cluster{si}(1:end));
        plot(sig_cluster{si}(1:end)-300, ones(1,segLength).*(yLimits(2)*0.5), 'LineWidth',2, 'Color', color{1}/255);
        hold on % plot multiple clusters
    end

    % save: bootstrapStatistics, lowerBound, upperBound, rate, rateSummary.
    % , rtSummary, rtTrialwise
    save(strcat(summaryName, '.mat'), 'rate', 'rateSummary', 'rtSummary', 'rtTrialwise', 'rtUnfiltered', 'trialData');

    % Set the PaperOrientation property to landscape
    set(gcf, 'PaperOrientation', 'landscape');

    % Get the current figure position in pixels
    fig_position_pixels = f.Position;
    
    % Get the screen resolution in pixels per inch
    screen_ppi = get(0, 'ScreenPixelsPerInch');
    
    % Convert figure position from pixels to inches
    fig_position_inches = fig_position_pixels(3:4) / screen_ppi;
    
    % Set the PaperPosition property to fit the entire figure
    set(f, 'PaperPosition', [0 0 fig_position_inches/2]);
        
    saveas(gcf, strcat(figName, '.fig')); % Save as FIG format
    saveas(gcf, strcat(figName, '.pdf')); % Save as PDF format
    %print(strcat(figName, '.pdf'), '-dpdf', '-bestfit'); % works but lines
    %are too thick

    close all
%     %% histogram of all RTs
%     figure
%     histogram(RTcard_all ,'FaceColor', color{1}/255, 'FaceAlpha', 0.5)
%     hold on
%     histogram(RTobl_all, 'FaceColor', color{2}/255, 'FaceAlpha', 0.5)
%     title(sprintf('RT Range for %s', subjects{ii}))
%     percentile_95 = prctile([RTcard_all; RTobl_all], 95);
end


%% Computing confidence / significance

function [lowerBound, upperBound] = findCI(bootstrapStatistics, percentileRange)

    if percentileRange == 68
        lPerc = 16; uPerc = 84;
    elseif percentileRange == 95
        lPerc = 2.5; uPerc = 97.5;
    end

    % Calculate 68% confidence intervals for each time point
    lowerBound = prctile(bootstrapStatistics, lPerc, 1);
    upperBound = prctile(bootstrapStatistics, uPerc, 1);

end

function [bootstrapStatistics] = bootstrapData(rate, fieldName)

    % Number of bootstrap samples
    numBootstrapSamples = 1000;
    
    % Preallocate array for storing statistics for each bootstrap sample
    bootstrapStatistics = zeros(numBootstrapSamples, size(rate.(fieldName), 2));
    
    % Bootstrap resampling
    for i = 1:numBootstrapSamples
        % Generate a bootstrap sample by sampling rows with replacement
        bootstrapSampleIndices = randi(size(rate.(fieldName), 1), size(rate.(fieldName), 1), 1);
        bootstrapSample = rate.(fieldName)(bootstrapSampleIndices, :);
        
        % Compute the statistic of interest for the bootstrap sample
        % For example, you can use mean, median, etc.
        bootstrapStatistic = nanmean(bootstrapSample, 1)*60; % Compute mean for each time point
        
        % Store the statistic for this bootstrap sample
        bootstrapStatistics(i, :) = bootstrapStatistic;
    end

end

function sig_cluster = checksignificance_perm(rateSummary)
        windowSize = 10; 
        
        % for permutation testing
        dependent_samples = false; % false because we are not directly comparing sessions (sometimes different # of sessions) 
        p_threshold = 0.05;
        num_permutations = 1000; 
        two_sided = false;
        num_clusters =[];

        fields = fieldnames(rateSummary);
        fields{1}
        fields{2}
        data1 = rateSummary.(fields{1});
        data2 = rateSummary.(fields{2});
    
        % Apply the moving average filter to each row separately
        smoothedData = zeros(size(data1));
        smoothedData2 = zeros(size(data2));

        for i = 1:size(data1, 1)
            smoothedData1(i, :) = movmean(data1(i, :), windowSize);
        end
        for i = 1:size(data2, 1)
            smoothedData2(i, :) = movmean(data2(i, :), windowSize);
        end
    
        [clusters, p_values, ~, ~ ] = permutest( smoothedData1', smoothedData2', dependent_samples, ...
        p_threshold, num_permutations, two_sided, num_clusters );
        sig_cluster = clusters(p_values<p_threshold); % find significant cluster(s)

end

