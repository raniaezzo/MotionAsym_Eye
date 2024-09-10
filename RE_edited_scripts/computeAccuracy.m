clc; clear all; close all;
rng(0);

cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
load('uniqueCombinations.mat')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 

totaln = nan(40,length(subjects));
acc_count_ms = nan(40,length(subjects));
acc_count_noms = nan(40,length(subjects));


% check whether a MS occurs within a critical window of 1300-1800
% Define the range of columns to check
startCol = 1300;
endCol = 1800;

%%

for sub=1:length(subjects)
    
% Define the subject ID as a string
    subjectID = subjects{sub};

    % Define the main directories
    baseDir = [datadir, subjectID, '/ProcessedData/'];
    folders = {'Full_distance_radialtangential', 'Full_distance_non_radialtangential'};

    % Initialize cell arrays to store the concatenated matrices
    concatenatedMicrosaccadeMatrices = {};
    concatenatedEventMatrices = {};

    % Loop through the main folders
    for i = 1:length(folders)
        % Define the eyedata/microsaccades subfolder
        subFolder = [baseDir, folders{i}, '/eyedata/microsaccades/'];
        subFolder2 = [baseDir, folders{i}, '/eyedata/MATs/'];

        % Get the list of all files in the subfolder
        allFiles = dir(subFolder2);

        % Extract unique identifiers (e.g., 'HL') by finding 'XX######_tab.mat' files
        uniqueIdentifiers = {};
        for j = 1:length(allFiles)
            [~, fileName, ~] = fileparts(allFiles(j).name);
            if contains(fileName, '_tab')
                uniqueIdentifiers{end+1} = fileName(1:8); % Extract the first 2 characters as the identifier
            end
        end

        if ~isempty(uniqueIdentifiers)
            % Remove duplicates
            uniqueIdentifiers = unique(uniqueIdentifiers);

            % Loop through each unique identifier
            for k = 1:length(uniqueIdentifiers)
                identifier = uniqueIdentifiers{k};

                % Load the corresponding events file
                eventsFilePattern = fullfile(subFolder, ['*', identifier(1:2), '_events.mat']);
                eventsFile = dir(eventsFilePattern);
                if ~isempty(eventsFile)
                    eventsData = load(fullfile(subFolder, eventsFile(1).name));
                    eventsMatrix = eventsData.('EVENTS'); % Load the matrix


                end

                % Load the microsaccadeMatrix file
                microsaccadeFile = fullfile(subFolder2, [identifier, '_tab.mat']);
                if exist(microsaccadeFile, 'file')
                    microsaccadeData = load(microsaccadeFile);
                    microsaccadeMatrix = microsaccadeData.('tab'); % Load the matrix

                end

                eventsSize = size(eventsMatrix,1);
                tabSize = size(microsaccadeMatrix,1);

                cutSize = min(eventsSize, tabSize);

                %if ~isempty(k) && ~isempty(j)
                    % Store the events matrix
                    concatenatedEventMatrices{end+1} = eventsMatrix(1:cutSize,:);

                    % Store the microsaccade matrix
                    concatenatedMicrosaccadeMatrices{end+1} = microsaccadeMatrix(1:cutSize,:);
                %end

            end
        end
    end

    % Optional: Concatenate all loaded matrices (for example, vertically)
    if ~isempty(concatenatedMicrosaccadeMatrices)
        finalMicrosaccadeMatrix = vertcat(concatenatedMicrosaccadeMatrices{:});
    end

    if ~isempty(concatenatedEventMatrices)
        finalEventMatrix = vertcat(concatenatedEventMatrices{:});
    end

    %%
    % Check each row for at least one '1' in the specified columns
    hasOne = any(finalEventMatrix(:, startCol:endCol), 2);
    tokenStimulus = finalMicrosaccadeMatrix(:,10:11);

    % Find unique combinations of values across the three columns
    %[uniqueCombinations, ~, idx] = unique(tokenStimulus, 'rows');

    % Initialize matrixSelect to store logical columns for each unique combination
    matrixSelect = false(size(tokenStimulus, 1), size(uniqueCombinations, 1));

    % Initialize an index vector to store the matching indices
    idx = zeros(size(tokenStimulus, 1), 1);

    % Loop through each row of tokenStimulus and find the matching index in uniqueCombinations
    for i = 1:size(tokenStimulus, 1)
        % Find the row in uniqueCombinations that matches the current row of tokenStimulus
        match = ismember(uniqueCombinations, tokenStimulus(i, :), 'rows');

        % Store the index of the matching row (if found)
        if any(match)
            idx(i) = find(match);
        else
            idx(i) = NaN; % If no match is found, you can assign NaN or any other value
        end
    end
    
    % Create logical columns for each unique combination
    for i = 1:size(uniqueCombinations, 1)
        matrixSelect(:, i) = idx == i;
    end

    for si=1:40
       msCurrent = matrixSelect(:,si) & hasOne;
       nomsCurrent = matrixSelect(:,si) & ~hasOne; 

       % Count the number of 1s in each vector
        numOnesMs = sum(msCurrent);
        numOnesNoms = sum(nomsCurrent);

        % Determine the target number of 1s (the smaller of the two counts)
        targetOnes = min(numOnesMs, numOnesNoms);

        % Adjust msCurrent if it has more 1s
        if numOnesMs > targetOnes
            % Find the indices of the 1s in msCurrent
            msIndices = find(msCurrent);

            % Randomly select indices to change from 1 to 0
            indicesToChange = randsample(msIndices, numOnesMs - targetOnes);

            % Change the selected 1s to 0s
            msCurrent(indicesToChange) = 0;
        end

        % Adjust nomsCurrent if it has more 1s
        if numOnesNoms > targetOnes
            % Find the indices of the 1s in nomsCurrent
            nomsIndices = find(nomsCurrent);

            % Randomly select indices to change from 1 to 0
            indicesToChange = randsample(nomsIndices, numOnesNoms - targetOnes);

            % Change the selected 1s to 0s
            nomsCurrent(indicesToChange) = 0;
        end

        totaln(si,sub) = targetOnes;
        acc_count_ms(si,sub) = sum(finalMicrosaccadeMatrix(msCurrent, 14));
        acc_count_noms(si,sub) = sum(finalMicrosaccadeMatrix(nomsCurrent, 14));
    end
end


%%

% Define the sets of values for the conditions
dirCardinal = [5, 6, 7, 8];
dirOblique = [1, 2, 3, 4];
tiltLarge = [8];
tiltSmall = [0.5, 1, 2, 4];

% Create a boolean column vector based on the conditions
booleanCardLarge = ismember(uniqueCombinations(:, 1), dirCardinal) & ...
                ismember(uniqueCombinations(:, 2), tiltLarge);
            
booleanCardSmall = ismember(uniqueCombinations(:, 1), dirCardinal) & ...
                ismember(uniqueCombinations(:, 2), tiltSmall);
            
booleanOblLarge = ismember(uniqueCombinations(:, 1), dirOblique) & ...
                ismember(uniqueCombinations(:, 2), tiltLarge);

booleanOblSmall = ismember(uniqueCombinations(:, 1), dirOblique) & ...
                ismember(uniqueCombinations(:, 2), tiltSmall);

            
            %%
            
n = length(subjects);
ones_vector = zeros(n, 1);
jittered_vector = ones_vector + (rand(n, 1) - 0.5) * 0.25;
            
% return vector for unique combinations that are cardinal vs. oblique
wMScard_small = sum(acc_count_ms(booleanCardSmall,:),1)./sum(totaln(booleanCardSmall,:),1);
wMScard_large = sum(acc_count_ms(booleanCardLarge,:),1)./sum(totaln(booleanCardLarge,:),1);
wMSobli_small = sum(acc_count_ms(booleanOblSmall,:),1)./sum(totaln(booleanOblSmall,:),1);
wMSobli_large = sum(acc_count_ms(booleanOblLarge,:),1)./sum(totaln(booleanOblLarge,:),1);
 
woMScard_small = sum(acc_count_noms(booleanCardSmall,:),1)./sum(totaln(booleanCardSmall,:),1);
woMScard_large = sum(acc_count_noms(booleanCardLarge,:),1)./sum(totaln(booleanCardLarge,:),1);
woMSobli_small = sum(acc_count_noms(booleanOblSmall,:),1)./sum(totaln(booleanOblSmall,:),1);
woMSobli_large = sum(acc_count_noms(booleanOblLarge,:),1)./sum(totaln(booleanOblLarge,:),1);

% wMS = sum(acc_count_ms,2)./sum(totaln,2);
% woMS = sum(acc_count_noms,2)./sum(totaln,2);

% combine all cardinal, and all oblique
color = {[17, 119, 51],[51, 34, 136]}; 

indvSize = 450;
meanSize = 550;

figure
hold on
plot([1 2], [mean(wMScard_large), mean(woMScard_large)], 'k', 'LineWidth', 2)
scatter(ones(length(subjects),1)+jittered_vector, wMScard_large, indvSize, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(1, mean(wMScard_large), meanSize,'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(2*ones(length(subjects),1)+jittered_vector, woMScard_large, indvSize, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(2, mean(woMScard_large), meanSize,'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
plot([3 4], [mean(wMScard_small), mean(woMScard_small)], 'k', 'LineWidth', 2)
hold on
scatter(3*ones(length(subjects),1)+jittered_vector, wMScard_small, indvSize, 'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(3, mean(wMScard_small), meanSize,'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(4*ones(length(subjects),1)+jittered_vector, woMScard_small, indvSize, 'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(4, mean(woMScard_small), meanSize,'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
plot([5 6], [mean(wMSobli_large), mean(woMSobli_large)], 'k', 'LineWidth', 2)
hold on
scatter(5*ones(length(subjects),1)+jittered_vector, wMSobli_large, indvSize, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(5, mean(wMSobli_large), meanSize,'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(6*ones(length(subjects),1)+jittered_vector, woMSobli_large, indvSize, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(6, mean(woMSobli_large), meanSize,'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
plot([7 8], [mean(wMSobli_small), mean(woMSobli_small)], 'k', 'LineWidth', 2)
hold on
scatter(7*ones(length(subjects),1)+jittered_vector, wMSobli_small, indvSize, 'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(7, mean(wMSobli_small), meanSize,'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(8*ones(length(subjects),1)+jittered_vector, woMSobli_small, indvSize, 'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(8, mean(woMSobli_small), meanSize,'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on

% % Plot the identity line (y = x)
% minVal = min([wMS, woMS]); % Determine the minimum value for the line
% maxVal = max([wMS, woMS]); % Determine the maximum value for the line
% 
% plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.5);

xlim([0 9])
ylim([0.5 1.15])
ylabel('% correct')

yticks([0.5 0.6 0.7 0.8 0.9 1])
xticks(1:8)
xticklabels({'w', 'w/o', 'w', 'w/o', 'w', 'w/o', 'w', 'w/o'})

set(gca, 'LineWidth', 2);

f1 = gcf;
f1.Position = [334 778 889 513];
dpi = get(0, 'ScreenPixelsPerInch');
set(gca, 'FontName', 'Arial', 'FontSize', (513/dpi)*4.3);

pathtoSave = '/Users/rje257/Desktop/oculomotor_manuscript_figures/';
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gcf, 'PaperPositionMode', 'auto'); 

set(gca, 'LineWidth', 1.5);
print(gcf, '-dpdf', fullfile(pathtoSave,'fig8.pdf'));  % for PDF
