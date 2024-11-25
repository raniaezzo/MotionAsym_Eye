% create simulated MS data
clear all;
%close all;


MSmatrix = zeros(300, 2500);

t = linspace(0, 1000, 1000);
suppressionStart = 250;
reboundStart = 850;
stimOn = 250;

%signalOut = createRateSignal(suppressionStart, reboundStart, stimOn)
signalOut = createRateSignal2(suppressionStart, reboundStart, stimOn,t)

% Plot the original and smoothed signal
plot(t, signalOut, 'r', 'LineWidth',5);
hold on
xline(250, 'LineWidth',5, 'Color', [0.8 0.8 0.8]);
hold on
%scatter(reboundStart, signalOut(reboundStart), 450, 'MarkerFaceColor', [0.7 0.7 0], 'MarkerEdgeColor', 'w', 'LineWidth', 2);
%hold on
xline(750, 'LineWidth',5, 'Color', [0.8 0.8 0.8]);
%xlabel('Time (ms)');
%ylabel('Signal Amplitude');
xticks([])
yticks([])
f1 = gcf;
f1.Position = [1000 1042 957 295];
ylim([0 2])
box off

%%


% Adjust the signal to ensure max value corresponds to 1.5% probability
signal = signalOut / 2 * 0.002;  % Scale to ensure max of signal corresponds to 1.5% probability

% Create a raster matrix with 500 trials and 1000 ms
num_trials = 5000;
raster_matrix = zeros(num_trials, 1000);

% Generate random values and compare with the adjusted signal to get 0s and 1s
for i = 1:num_trials
    random_vals = rand(1, 1000);
    raster_matrix(i, :) = random_vals < signal;
end

% add random 1s
% Assume raster_matrix is your existing binary matrix (500 trials, 1000 ms)
[num_trials, num_timepoints] = size(raster_matrix);

% Define how many random 1s you want to add
num_random_ones = num_trials/10;  % You can change this as needed

% Randomly select positions to insert 1s
random_rows = randi([1 num_trials], 1, num_random_ones);  % Random row indices
random_cols = randi([1 num_timepoints], 1, num_random_ones);  % Random column indices

% Set the randomly selected positions to 1
for k = 1:num_random_ones
    raster_matrix(random_rows(k), random_cols(k)) = 1;
end
% done

% Pad each 1 with a random number of 1s (uniformly distributed between 50 and 300) to the immediate left and right row-wise, independently
padded_matrix = raster_matrix;  % Create a copy to modify

for i = 1:num_trials
    % Find indices of all '1's in the original raster_matrix (before modification)
    ones_idx = find(raster_matrix(i, :) == 1);
    
    for j = 1:length(ones_idx)
        current_idx = ones_idx(j);  % Get the index of the original 1 (TAG)
        
        % Determine start and end indices for padding (30 columns before and after)
        start_idx = max(1, current_idx - 15);  % Ensure start does not go below 1
        end_idx = min(1000, current_idx + 15);  % Ensure end does not exceed 1000
        
        % Set all columns in this range to 1
        padded_matrix(i, start_idx:end_idx) = 1;
    end
end

figure
imagesc(padded_matrix)


%%

% create movie

[nrows, ncols] = size(padded_matrix);

figure;
ylim([0 nrows])
xlim([0 ncols])
axis square

% Number of frames for the video
num_frames = nrows;

% Initialize matrix for falling simulation
falling_matrix = padded_matrix;

% MP4 settings
filename = 'falling_onesNEW2.mp4';
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = 125;  % Adjust the frame rate (you can change this as needed)
open(v);

for t = 1:num_frames
    % For each row (starting from the second to last row), check if "1"s can fall
    for row = size(falling_matrix, 1)-1:-1:1  % From second-to-last row upward
        for col = 1:size(falling_matrix, 2)
            if falling_matrix(row, col) == 1 && falling_matrix(row+1, col) == 0
                % Swap positions (1 falls down)
                falling_matrix(row, col) = 0;
                falling_matrix(row+1, col) = 1;
            end
        end
    end
    
    % Display the updated matrix
    imagesc(falling_matrix);
    %colormap(gray);  % Set colormap to grayscale for better visibility
    colormap(flipud(gray));
    axis off;
    ylim([1, 5500])
    %xlim([500,2250])
    title(['Frame ' num2str(t)]);
    drawnow;
    
    % Capture the frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);  % Write the frame to the video
end

% Close the video file
close(v);

% Display final "pile" in histogram
pile_heights = sum(falling_matrix == 1, 1);
figure;
bar(pile_heights);
%xlabel('Column');
%ylabel('Pile Height');
%title('Final Pile Height Histogram');



% Load the original video
videoObj = VideoReader(filename);  % Replace with your video file

% Create a VideoWriter object for the new, faster video
outputVideo = VideoWriter('falling_onesSPEDUP.mp4', 'MPEG-4');
outputVideo.FrameRate = videoObj.FrameRate;  % Keep original frame rate

% Open the VideoWriter object to write frames
open(outputVideo);

% Calculate the speed-up factor (120 seconds -> 4 seconds means 30x speed-up)
speedupFactor = 30;

% Read frames from the input video and write the sped-up frames to the new video
frameCount = 0;
while hasFrame(videoObj)
    frame = readFrame(videoObj);
    frameCount = frameCount + 1;
    
    % Write every Nth frame to speed up the video by the speedup factor
    if mod(frameCount, speedupFactor) == 0
        writeVideo(outputVideo, frame);
    end
end

% Close the VideoWriter object
close(outputVideo);

disp('Video has been sped up and saved!');








%%


function signalOut = createRateSignal2(suppressionStart, reboundStart, stimOn, t)
    
if suppressionStart > stimOn
    error("suppression must occur at or before stimOn")
end

% Initialize the signal with ones
signalOut = ones(size(t));

trialLength = length(t);

% hovering baseline
initialSegTime = 1:suppressionStart;

if length(initialSegTime) == 1
    % this should be fixed (only shifted)
    initialSegValues = [];
else
    initialSegValues = ones(1,length(initialSegTime));
end

% how long is the suppression ramp
minSupDuration = 200;
% what time does the suppression stop decreasing
suppDone = suppressionStart+minSupDuration+(suppressionStart-stimOn);
suppressionrampTime = suppressionStart:suppDone;
%suppressionrampValues = linspace(1,0, length(suppressionrampTime)-1); 
suppressionrampTime = linspace(0, pi, length(suppressionrampTime)-1);  % Time over pi for cosine
suppressionrampValues = 0.5 * (1 + cos(suppressionrampTime));  % Cosine ramp from 1 to 0

% how long is the signal hovering 0
minTime = suppDone:reboundStart;

if length(minTime) == 1
    % this should be fixed (only shifted)
    minTimeValues = [];
else
    minTimeValues = zeros(1,length(minTime)-1);
end

returnBaseline = reboundStart:trialLength;

ramp_duration = length(returnBaseline)-1; %300;
ramp_time = linspace(0, 2*pi*.75, ramp_duration);  % Create a time vector for the cosine ramp
ramp_timeVals = 1 + (cos(ramp_time) * -(2 - 1));

signalOut = [initialSegValues, suppressionrampValues, minTimeValues, ramp_timeVals]; 

%     Fs = 1000;  % Sampling frequency
%     t = linspace(0, 1, Fs);  % Time vector for 1 second
% 
%     % Generate random noise
%     noise = randn(size(t));
% 
%     % Design a low-pass filter to create low-frequency noise
%     Fc = 3;  % Cutoff frequency for low-frequency noise (5 Hz)
%     [b, a] = butter(2, Fc/(Fs/2), 'low');  % 2nd-order Butterworth low-pass filter
%     low_freq_noise = filtfilt(b, a, noise);  % Filter the random noise
% 
%     % Add the low-frequency noise to the signal
%     signalOut = signalOut + 0.2 * low_freq_noise;  % Scale noise as needed
%     signalOut = abs(signalOut);

end






















function signalOut = createRateSignal(suppressionStart, reboundStart, stimOn)

    Fs = 1000;  % Sampling frequency in Hz
    Fc = 8;    % Cutoff frequency (adjust as needed)

    if suppressionStart == stimOn
        suppressionStart = suppressionStart+(Fs/Fc); % adjust to account for low pass
    end

    figure
    % Create a time vector from 0 to 1000 ms
    t = linspace(0, 1000, 1000);
    
    % Initialize the signal with ones
    signalOut = ones(size(t));
    
    % Smoothly vary around 1 before suppressionStart
    signalOut(1:suppressionStart) = 1 + 0.1 * sin(2 * pi * 0.01 * t(1:suppressionStart));
    
    % Suppress signal from suppressionStart to stimOn
    signalOut(suppressionStart:stimOn) = linspace(1, 0, stimOn - suppressionStart + 1);
    
    % After stimOn, vary around 0 (but never go below)
    signalOut(stimOn:reboundStart) = 0.1 * abs(sin(2 * pi * 0.01 * t(stimOn:reboundStart)));
    
    % From reboundStart, signal increases gradually to peak at 2
    signalOut(reboundStart:end) = linspace(0, 2, length(t(reboundStart:end)));
    
    % Smoothly return to 1 before the end of the signal
    signalOut(end-100:end) = linspace(2, 1, 101);

    alpha = 0.1;  % Smoothing factor (0 < alpha <= 1)
    initial_value = 1;
    signalOut = filter(alpha, [1, -(1 - alpha)], signalOut, initial_value * (1 - alpha));

    % From reboundStart, signal follows a cosine ramp that peaks at 2 and returns to 1
    ramp_duration = length(t(reboundStart:end));
    ramp_time = linspace(0, pi, ramp_duration);  % Create a time vector for the cosine ramp
    
    signalOut(reboundStart:end) = 1 + (cos(ramp_time) * -(2 - 1));  % Cosine ramp, peaks at 2 and returns to 1

%     % Define Gaussian kernel
%     sigma = 20;  % Standard deviation
%     gaussian_kernel = fspecial('gaussian', [1, 100], sigma);
% 
%     % Apply padding (reflection)
%     padding_size = 50;  % Adjust based on kernel size
%     padded_signal = [flip(signalOut(1:padding_size)), signalOut, flip(signalOut(end-padding_size+1:end))];
% 
%     % Convolve the padded signal with the Gaussian kernel
%     smoothed_signal = conv(padded_signal, gaussian_kernel, 'same');
% 
%     % Remove padding
%     signalOut = smoothed_signal(padding_size+1:end-padding_size);


end

