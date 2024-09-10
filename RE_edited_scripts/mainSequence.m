% for pupil, retrieve values for current trial and next trial from ASC.
%clc; clear all; close all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 

%%
% Define the subject ID as a string
subjectID = 'S01';

% Define the directories
dir1 = fullfile(datadir, 'Data_DI_wEYE', subjectID, 'ProcessedData', 'Full_distance_radialtangential', 'eyedata', 'microsaccades');
dir2 = fullfile(datadir, 'Data_DI_wEYE', subjectID, 'ProcessedData', 'Full_distance_non_radialtangential', 'eyedata', 'microsaccades');

% Get list of files ending in _summary.mat in both directories
filesDir1 = dir(fullfile(dir1, '*microsaccadeMatrix.mat'));
filesDir2 = dir(fullfile(dir2, '*microsaccadeMatrix.mat'));

% Initialize an empty cell array to store concatenated matrices
concatenatedMatrices = {};

% Load and concatenate matrices from the first directory
for i = 1:length(filesDir1)
    % Load the .mat file
    data = load(fullfile(dir1, filesDir1(i).name));
    
    % Assume the variable name is consistent and known, e.g., 'dataMatrix'
    varName = fieldnames(data);
    matrix = data.(varName{1}); % Load the matrix
    
    % Assign the matrix to a new variable name to avoid overwriting
    concatenatedMatrices{end+1} = matrix;
end

% Load and concatenate matrices from the second directory
for i = 1:length(filesDir2)
    % Load the .mat file
    data = load(fullfile(dir2, filesDir2(i).name));
    
    % Assume the variable name is consistent and known, e.g., 'dataMatrix'
    varName = fieldnames(data);
    matrix = data.(varName{1}); % Load the matrix
    
    % Assign the matrix to a new variable name to avoid overwriting
    concatenatedMatrices{end+1} = matrix;
end

% Convert the cell array of matrices into a single concatenated matrix
% Assuming horizontal concatenation; adjust as needed
finalMatrix = vertcat(concatenatedMatrices{:});





%%


figure
scatter(finalMatrix(:,11),finalMatrix(:,3), 450, 'o', 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none','MarkerFaceAlpha', 0.05) %0.008)
ylabel('peak velocity (deg/s)')
xlabel('amplitude (deg)')

set(gca, 'FontName', 'Arial', 'FontSize', 20); % this changes
set(gca, 'LineWidth', 4);

f1 = gcf;
f1.Position = [82 874 397 285];

xticks([0 .25 .5 .75 1])

% CANNOT SAVE IN VECTOR FORMAT--TOO MANY POINTS!
% to adjust ratio / font
% Get current figure properties
% Get current figure size (width in inches)
current_size = get(gcf, 'PaperPosition');  % [left, bottom, width, height]
width = current_size(3);                  % Width of the figure
font_size = get(gca, 'FontSize');
current_ratio = font_size / width;

screen_size = get(0, 'ScreenSize');  % Get screen size [left, bottom, width, height]
set(gcf, 'Position', [1 1 screen_size(3)/2, (screen_size(3)/2)*(2/3)]); %screen_size/2);   % Set figure position to fill the screen

% Adjust font size based on the new width
dpi = get(0, 'ScreenPixelsPerInch');
screen_size_inches = (screen_size/2) / dpi;
new_font_size = current_ratio * screen_size_inches(3);
set(gca, 'FontSize', new_font_size);


