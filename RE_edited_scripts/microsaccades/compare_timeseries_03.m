clc; clear all; %close all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 

analysis_type = 'tilt'; %'direction'; 
trialwise = 1;
percentileRange = 68; % for plotting, do 68% CIs
counter = 1; % for figures

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

%% MS latency: Probability density (z-score)

% calculate the z-score from the timepoints of the first before suppression
% to first MS after suppression

stimonset = 300; % this data has 1000 ms cutoff from beginning
stimoffset = 800; 

for ii = 1:1

    figure
    summaryMSPath = fullfile(datadir, subjects{ii}, 'ProcessedData', 'Summary', 'microsaccades');
    fileName = sprintf('%s_%s_%s_%s_%i_%i', subjects{ii}, analysis_type, fieldNames{1}, fieldNames{2}, trialwise, percentileRange);
    summaryName = fullfile(summaryMSPath, fileName);
    load(strcat(summaryName, '.mat'))

    % find the suppression point (during stimulus)
    combinedData = nanmean([rate.(fieldNames{1}); rate.(fieldNames{2})]);
    [~,suppressionTime] = min(combinedData(stimonset:stimoffset)); % stim on to offset (adds back offset in the function)

    % this is for all trials (TO DO: separate by accuracy)
    %[postOnsetMS_cond1, latencyCond1, RT_outputCond1, accuracy_Cond1] = compute_latency(rate, '', suppressionTime, rtTrialwise, trialData);

    figure(counter);

    % all data
    [postOnsetMS_cond, latencyCond, RT_outputCond, accuracy_Cond, selectedTrials] = compute_latency(rate, '', suppressionTime, rtTrialwise, trialData, stimonset);

    subplot(1,2,1)
    postOnsetMS_corr = postOnsetMS_cond(accuracy_Cond == 1);
    postOnsetMS_incorr = postOnsetMS_cond(accuracy_Cond == 0);
    RT_outputCond_corr = RT_outputCond(accuracy_Cond == 1);
    RT_outputCond_incorr = RT_outputCond(accuracy_Cond == 0);

    plot(RT_outputCond_corr, postOnsetMS_corr, '.', 'Color', 'g')
    hold on
    [x_fit, y_fit] = fitLine(RT_outputCond_corr, postOnsetMS_corr);
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'g');
    hold on
    [x_fit, y_fit] = fitLine(RT_outputCond_incorr, postOnsetMS_incorr);
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'r');
    plot(RT_outputCond_incorr, postOnsetMS_incorr, '.', 'Color', 'r')
    hold on

    correlation_coefficient1 = corr(RT_outputCond_corr, postOnsetMS_corr', 'type', 'Spearman', 'rows', 'complete')
    correlation_coefficient2 = corr(RT_outputCond_incorr, postOnsetMS_incorr', 'type', 'Spearman', 'rows', 'complete')

    xlim([500 2000])
    xlabel('reaction time (ms)')
    ylabel('onset of first microsaccade (ms)')

    % global stats to plot probablity density function
    avTotal = mean(latencyCond);
    stdevTotal = std(latencyCond);

    subplot(1,2,1)
    counter = counter+1;

    for fi=1:length(fieldnames(rate))

        [postOnsetMS_cond, latencyCond, RT_outputCond, accuracy_Cond, qualifyingTrials, dataCheck] = compute_latency(rate, fieldNames{fi}, suppressionTime, rtTrialwise, trialData, stimonset);
        
        figure(counter-1);
        % plot(RT_outputCond1, latencyCond1, '.b')
        subplot(1,2,2)
        plot(RT_outputCond, postOnsetMS_cond, '.', 'Color', color{fi}/255)
        hold on
        xlim([500 2000])
        xlabel('reaction time (ms)')
        ylabel('onset of first microsaccade (ms)')
        %ylabel('latency (ms)')
        % correlation_coefficient1 = corr(RT_outputCond1, latencyCond1', 'type', 'Spearman')
        correlation_coefficient1 = corr(RT_outputCond, postOnsetMS_cond', 'type', 'Spearman', 'rows', 'complete')
        
        [x_fit, y_fit] = fitLine(RT_outputCond, postOnsetMS_cond);
        plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', color{fi}/255);

        % Density plot (need to fix
        figure(counter);
        z_scored_data = (latencyCond - avTotal) / stdevTotal;
        x_values = linspace(-5, 5, 100);
        [f, xi] = ksdensity(latencyCond, x_values);
        [fnorm, xinorm] = ksdensity(z_scored_data, x_values); % Calculate kernel density estimate


        subplot(2, 2, 1)
        histogram(latencyCond, 50, 'FaceColor', color{fi}/255, 'EdgeColor', 'none')
        hold on
        xline(median(latencyCond), 'Color', color{fi}/255, 'LineWidth', 2)
        subplot(2, 2, 2)
        plot(xi, f, 'Color', color{fi}/255); % Plot the KDE
        hold on
        subplot(2, 2, 3)
        plot(xinorm, fnorm, 'Color', color{fi}/255); % Plot the KDE
        hold on
        xline(median(z_scored_data), 'Color', color{fi}/255)
        xlabel('latency (z-score)');
        ylabel('probability density');
        title('Kernel Density Estimate');
        subplot(2, 2, 4)
        h = cdfplot(latencyCond);
        set(h, 'Color', color{fi}/255)
        hold on
    end

    counter = counter+1;



end



% av = mean([latencyCond, latencyCond2]);
% stdev = std([latencyCond, latencyCond2]);
% z_scored_data1 = (latencyCond - av) / stdev;
% z_scored_data2 = (latencyCond2 - av) / stdev;

% kernel density
% Define the range of points
% x_values = linspace(-5, 5, 100);
% [f1, xi1] = ksdensity(z_scored_data1, x_values); % Calculate kernel density estimate
% [f2, xi2] = ksdensity(z_scored_data2, x_values);

% do it again it time units
% x_values = linspace(0, 1600, 100);
% [f3, xi3] = ksdensity(latencyCond, x_values); % Calculate kernel density estimate
% [f4, xi4] = ksdensity(latencyCond2, x_values);

% figure
% subplot(2, 2, 1)
% histogram(latencyCond, 50, 'FaceColor', 'blue', 'EdgeColor', 'none')
% hold on
% histogram(latencyCond2, 50, 'FaceColor', 'red', 'EdgeColor', 'none')
% hold on
% xline(median(latencyCond), 'b', 'LineWidth', 2)
% hold on
% xline(median(latencyCond2), 'r', 'LineWidth', 2)
% subplot(2, 2, 2)
% % Plot the KDE
% plot(xi3, f3, 'b');
% hold on
% plot(xi4, f4, 'r');
% subplot(2, 2, 3)
% plot(xi1, f1, 'b');
% hold on
% plot(xi2, f2, 'r');
% hold on
% xline(median(z_scored_data1), 'b')
% hold on
% xline(median(z_scored_data2), 'r')
% xlabel('latency (z-score)');
% ylabel('probability density');
% title('Kernel Density Estimate');
% subplot(2, 2, 4)
% cdfplot(latencyCond);
% hold on
% cdfplot(latencyCond2);

%% estimate the peak after inflection point

%% estimate the slope from inflection to peak

%% estimate the area under the curve from inflection to peak

%% estimate the area under the curve during stimulus period

%% estimate the shift at peak correation

%% estimate the horizontal shift (in time) from inflection point to peak

%% estimate the vertical shift (in rate) from inflection point to peak

%% Trying to find inflection point

data1 = nanmean(rate.cardinal,1)*60;
data2 = nanmean(rate.oblique,1)*60;

% estimate the inflection point 

    % smooth curve
    windowSize = 150; 

    % Apply the moving average filter to each row separately
    smoothedData = zeros(size(data1));
    smoothedData2 = zeros(size(data2));

    for i = 1:size(data1, 1)
        smoothedData1(i, :) = movmean(data1(i, :), windowSize);
    end
    for i = 1:size(data2, 1)
        smoothedData2(i, :) = movmean(data2(i, :), windowSize);
    end

    % find 0 derivative
    % Compute the derivative of the time series
    derivative = diff(smoothedData1);

    figure
    %subplot(2,1,1)
    plot(smoothedData1, 'b')
    hold on
    %plot(smoothedData2, 'r')
    %hold on
    xline(300)
    hold on
    xline(800)
    hold on
%     subplot(2,1,2)
%     plot(derivative)
%     hold on
%     xline(300)
%     hold on
%     xline(800)
% 
%     hold on
    
    % Compute the moving average using a window of 100 ms
    moving_average = movmean(data1, windowSize);
    
    % Compute the derivative of the moving average
    derivative = diff(moving_average);
    
    % Find the index where the derivative is maximum
    [max_derivative, max_derivative_index] = max(derivative(300:1000));
    
    % Find the starting point of the largest increase
    starting_point = max_derivative_index; % to account for the shift

    xline(starting_point, ':')
    
    % Find the indices where the derivative changes from negative to positive
    positive_derivative_indices = find(derivative > 0);
    positive_derivative_indices(1:300) = NaN;
    
    consecutive_thresh = 50;


    % Initialize variables
    consecutive_count = 0;
    consecutive_start = 0;

    % Find the starting index of at least 50 consecutive positive derivatives
    for i = 1:length(positive_derivative_indices)
        if consecutive_count >= consecutive_thresh
            break;  % Exit the loop if we have found at least N consecutive positive derivatives
        end
        
        if consecutive_start == 0
            consecutive_start = positive_derivative_indices(i);
            consecutive_count = 1;
        elseif positive_derivative_indices(i) == consecutive_start + consecutive_count
            consecutive_count = consecutive_count + 1;
        else
            consecutive_start = positive_derivative_indices(i);
            consecutive_count = 1;
        end
    end
    
    % If we found at least N consecutive positive derivatives, retrieve the corresponding values
    if consecutive_count >= consecutive_thresh
        positive_derivative_values = smoothedData1(consecutive_start:consecutive_start+consecutive_count-1);
    else
        positive_derivative_values = [];
        disp('There are less than N consecutive positive derivatives.');
    end

    plot(positive_derivative_indices, 0*ones(length(positive_derivative_indices), 1), 'or')
    [~, timeVal] = find(~isnan(positive_derivative_indices));
    xline(positive_derivative_indices(timeVal(1)), ':')