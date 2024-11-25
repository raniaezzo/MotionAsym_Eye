clc;
clear all;
close all;
rng(0)


% Total duration and frames
duration = 10;  % seconds
fps = 30;       % frames per second
num_frames = duration * fps;  % total number of frames

% ********** Smooth Gaze Movements Using Brownian Motion **********
% Initialize gaze positions with zeros
gaze_x_smooth = zeros(1, num_frames);
gaze_y_smooth = zeros(1, num_frames);

% Define the radius of the eyeball (limit for the iris)
eye_radius = 0.5;

% Brownian motion parameters
step_size = 0.04;  % Adjust step size for smooth motion (can adjust this)

% Generate random steps for Brownian motion
for i = 2:num_frames
    % Random step in x and y directions (small increments)
    gaze_x_smooth(i) = gaze_x_smooth(i-1) + step_size * randn();
    gaze_y_smooth(i) = gaze_y_smooth(i-1) + step_size * randn();
    
    % Calculate the distance of the new gaze position from the center (0, 0)
    distance_from_center = sqrt(gaze_x_smooth(i)^2 + gaze_y_smooth(i)^2);
    
    % Check if the new position is outside the boundary of the eyeball
    if distance_from_center > eye_radius
        % Scale the position back so that it lies within the boundary
        scale_factor = eye_radius / distance_from_center;
        gaze_x_smooth(i) = gaze_x_smooth(i) * scale_factor;
        gaze_y_smooth(i) = gaze_y_smooth(i) * scale_factor;
    end
end

% ********** Center Gaze Positions Around 0 **********
% Subtract the mean from the gaze positions to ensure the mean is 0
gaze_x_smooth = gaze_x_smooth - mean(gaze_x_smooth);
gaze_y_smooth = gaze_y_smooth - mean(gaze_y_smooth);

% ********** Jerky Movements **********
% We'll introduce jerky movements every 1 to 2 seconds
jerky_amplitude = 0.6;  % Amplitude of jerky movements (can adjust this)
jerk_intervals = randi([fps, 2*fps], 1, floor(duration/2));  % Random intervals for jerks
jerk_times = cumsum(jerk_intervals);  % Cumulative times for when jerks happen
jerk_times(jerk_times > num_frames) = [];  % Remove any jerks beyond the duration

% Apply jerky movements at specified intervals
for i = 1:length(jerk_times)
    % Pick a random angle for the jerky jump (to randomize direction)
    random_angle = 2 * pi * rand();
    jerk_x = jerky_amplitude * cos(random_angle);  % Jerky jump in x direction
    jerk_y = jerky_amplitude * sin(random_angle);  % Jerky jump in y direction
    
    % Override the smooth gaze with the jerky movement at the specified time
    jerk_frame = jerk_times(i);
    new_gaze_x = gaze_x_smooth(jerk_frame:end) + jerk_x;
    new_gaze_y = gaze_y_smooth(jerk_frame:end) + jerk_y;
    
    % Reapply boundary constraints to jerky movements
    for j = jerk_frame:num_frames
        distance_from_center = sqrt(new_gaze_x(j-jerk_frame+1)^2 + new_gaze_y(j-jerk_frame+1)^2);
        if distance_from_center > eye_radius
            scale_factor = eye_radius / distance_from_center;
            new_gaze_x(j-jerk_frame+1) = new_gaze_x(j-jerk_frame+1) * scale_factor;
            new_gaze_y(j-jerk_frame+1) = new_gaze_y(j-jerk_frame+1) * scale_factor;
        end
    end
    
    gaze_x_smooth(jerk_frame:end) = new_gaze_x;
    gaze_y_smooth(jerk_frame:end) = new_gaze_y;
end

% ********** Plot on Cartesian Coordinates **********
figure;
plot(gaze_x_smooth, gaze_y_smooth, 'b-', 'LineWidth', 2);  % Plot the combined x and y trajectory
hold on;
scatter(gaze_x_smooth(jerk_times), gaze_y_smooth(jerk_times), 50, 'r', 'filled');  % Highlight jerky movements

% Add labels and title
xlabel('Gaze X Position');
ylabel('Gaze Y Position');
axis equal;  % Ensure equal scaling on both axes
grid on;

% Mark the start and end points for clarity
scatter(gaze_x_smooth(1), gaze_y_smooth(1), 100, 'g', 'filled');  % Start point (green)
scatter(gaze_x_smooth(end), gaze_y_smooth(end), 100, 'k', 'filled');  % End point (black)

%legend({'Gaze Path', 'Jerky Movements', 'Start', 'End'});


% ********** Generate Pupil Sizes **********
% Parameters for sinusoidal oscillation
pupil_mean = 1.5;     % Mean pupil size
pupil_amplitude = 0.3; % Amplitude of pupil size change
pupil_freq = 1;       % Frequency of pupil oscillation (1 cycle per second)

% Time vector for sinusoidal pupil size oscillation
t = linspace(0, duration, num_frames);

% Generate smooth sinusoidal oscillation for pupil size
pupil_sizes = pupil_mean + pupil_amplitude * sin(2 * pi * pupil_freq * t);

% ********** Add Random Toggling **********
% Toggle probability and amplitude
toggle_prob = 0.05;  % Probability of a random toggle at each frame (5% chance)
toggle_amplitude = 0.2;  % Amplitude of the random toggles

amplitude_variation = 0.2 * randn(1, num_frames);  % Add variability to amplitude
pupil_sizes = pupil_mean + (pupil_amplitude + amplitude_variation) .* sin(2 * pi * pupil_freq * t);

% Add Random Toggling
for i = 1:num_frames
    if rand() < toggle_prob
        toggle = toggle_amplitude * (2 * rand() - 1);  % Random value in [-toggle_amplitude, toggle_amplitude]
        pupil_sizes(i) = pupil_sizes(i) + toggle;
    end
end

pupil_sizes = 0.5.*(pupil_sizes+[pupil_sizes(100:end), pupil_sizes(1:99)]);

% ********** Plot Pupil Sizes Over Time **********
figure;
plot(t, pupil_sizes, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Pupil Size');
title('Pupil Size Oscillation with Random Toggling');
grid on;

% Mark the min and max points for clarity
[~, min_index] = min(pupil_sizes);
[~, max_index] = max(pupil_sizes);
hold on;
scatter(t(min_index), pupil_sizes(min_index), 100, 'b', 'filled');  % Mark min
scatter(t(max_index), pupil_sizes(max_index), 100, 'r', 'filled');  % Mark max
%legend({'Pupil Size', 'Min Size', 'Max Size'});



% Now you can use these gaze_x_smooth and gaze_y_smooth in your rotating_eye function.
gaze_positions = [gaze_x_smooth', gaze_y_smooth'];  
gaze_positions = gaze_positions.* .3;

pupil_sizes;

Fs = 1000;  % Sampling frequency (Hz)
Fc = 50;    % Cutoff frequency (Hz)

% Design a low-pass filter
[b, a] = butter(2, Fc/(Fs/2), 'low');  % 2nd-order low-pass Butterworth filter

% Apply the filter to your signal
pupil_sizes = filtfilt(b, a, pupil_sizes);
 

%% Prepare jerk trace

% Initialize logical array for jerk movements (true for jerk, false for smooth)
jerk_trace = false(1, num_frames);

% Mark the jerk frames
for i = 1:length(jerk_times)
    % Specify how long the jerk movement lasts (e.g., for a fixed duration)
    jerk_duration = fps / 6;  % Example: jerk lasts for half a second (adjust as needed)
    
    % Mark the frames affected by the jerk
    jerk_start = jerk_times(i);
    jerk_end = min(jerk_start + jerk_duration - 1, num_frames);  % Ensure it doesn't exceed the total number of frames
    jerk_trace(jerk_start:jerk_end) = true;
end


%%
% % Call the function to display two eyes


rotating_eyes(gaze_positions, pupil_sizes, eye_radius, jerk_trace);

function rotating_eyes(gaze_positions, pupil_sizes, eye_radius, jerk_indices)

    % Create VideoWriter objects for each figure
    video1 = VideoWriter('video1.mp4', 'MPEG-4');
    video2 = VideoWriter('video2.mp4', 'MPEG-4');
    video3 = VideoWriter('video3.mp4', 'MPEG-4');
    video4 = VideoWriter('video4.mp4', 'MPEG-4');
    
    % Set frame rate (optional, default is 30)
    video1.FrameRate = 15;
    video2.FrameRate = 15;
    video3.FrameRate = 15;
    video4.FrameRate = 15;
    
    % Open the video files
    open(video1);
    open(video2);
    open(video3);
    open(video4);


    figure1 = figure(1);
    figure2 = figure(2);
    figure3 = figure(3);
    figure4 = figure(4);

    background_img = imread('/Users/rje257/Desktop/generic_image.jpeg'); 

    % Convert to double for calculation, then back to uint8 after adjusting contrast
    background_img = im2double(background_img);
    
    % Scale pixel values towards the mean
    mean_val = mean(background_img(:));  % Compute the overall mean
    scaling_factor = 0.5;  % Lower this value to reduce contrast more (0 to 1)
    
    % Adjust contrast by scaling values towards the mean
    low_contrast_img = mean_val + scaling_factor * (background_img - mean_val);
    
    % Convert back to uint8
    background_img = im2uint8(low_contrast_img);
    
    % Get the size of the image
    [img_height, img_width, ~] = size(background_img);
    
    % Crop the image to the largest square possible (centered)
    if img_width > img_height
        % If the width is larger, crop the sides
        crop_x_start = floor((img_width - img_height) / 2) + 1;
        crop_x_end = crop_x_start + img_height - 1;
        background_img_cropped = background_img(:, crop_x_start:crop_x_end, :);
    else
        % If the height is larger or equal, no need to crop the height (or crop vertically similarly)
        background_img_cropped = background_img;
    end
    %background_img_cropped = flipud(background_img_cropped);

    % new mini image
    [img_height, img_width, ~] = size(background_img_cropped);
    new_size = round(img_width / 4);
    start_x = round((img_width - new_size) / 2);
    start_y = round((img_height - new_size) / 2);

    background_img_cropped_zoom = background_img_cropped(start_y:(start_y + new_size - 1), start_x:(start_x + new_size - 1), :);

    fps = 30; 
    % Total duration and frames
    num_frames = size(gaze_positions, 1);  % Total number of frames

    % Create the figure window for the eyes
    figure(figure1);
    %subplot(1,2,1);  % Subplot for eye movements
    set(figure1, 'Color', 'w')
    axis equal;
    axis([-2 3 -1.5 1.5]);
    %axis([-3 3 -1.5 1.5]);  % Set axis limits
    hold on;
    axis off;  % Hide the axis
    
    % Subplot for the real-time gaze trajectory
    %subplot(1,2,2);  % Subplot for gaze positions
    figure(figure2);
    hold on;
    xlabel('Gaze X Position', 'FontSize', 18);
    ylabel('Gaze Y Position', 'FontSize', 18);
    title('Gaze Positions with Jerky Movements');
    set(figure2, 'Color', 'w')
    axis equal;
    axis([-1 1 -1 1]);  % Adjust based on the gaze limits
    xticks([])
    yticks([])
    grid on;
    
    %Display the image in the background, scaled to fit [-1, 1] range
    imagesc([-1, 1], [-1, 1], background_img_cropped);
    hold on;  % Ensure plotting happens on top of the image
    
    % Set the axis limits to match the normalized range
    axis([-1, 1, -1, 1]);
    
    % Set the axis limits to match the centered image size
    set(gca,'YDir','reverse'); 

    % Eyeball (white circle)
    theta = linspace(0, 2*pi, 100);
    x_eye = cos(theta);
    y_eye = sin(theta);

    % Iris (brown oval)
    iris_radius_x = 0.5;
    iris_radius_y = 0.5;
    iris_color = [0.6 0.3 0];  % Brownish color

    % Default Pupil (black oval) size
    default_pupil_radius_x = 0.1;
    default_pupil_radius_y = 0.1;
    pupil_color = 'k';  % Black color

    figure(figure3);
    imagesc([-1, 1], [-1, 1], background_img_cropped);
    hold on;
    %title('Zoomed in');
    axis equal;
    axis([-.25 .25 -.25 .25]);  % Adjust based on the gaze limits
    set(figure3, 'Color', 'w')
    xticks([])
    yticks([])
    set(gca,'YDir','reverse'); 
    %grid on;

    figure(figure4)
    set(figure4, 'Color', 'w')
    ylim([1.1, 2])
    xlim([0 300])
    yticks([])
    xticks([])
    axis equal;

    % Initialize real-time plot for gaze positions
    figure(figure2);
    h_gaze = scatter(0, 0, 100, 'k+', 'LineWidth', 4);  % Blue for smooth trajectory

    %h_jerk = scatter(0, 0, 100, 'k+', 'LineWidth', 4); %scatter([], [], 50, 'r', 'filled');  % Red for jerky movements

    figure(figure3);
    h_gaze2 = scatter(0, 0, 400, 'k+', 'LineWidth', 4);

    % Loop through each gaze position and update the iris, pupil, and gaze trajectory
    for i = 1:num_frames
        % Get the current gaze position
        gaze_x = gaze_positions(i, 1);  % Gaze position on x-axis
        gaze_y = gaze_positions(i, 2);  % Gaze position on y-axis
        
        % Get the current pupil size
        pupil_scale = pupil_sizes(i);  % Scale factor for pupil size
        
        % Calculate the distance of the new gaze position from the center (0, 0)
        distance_from_center = sqrt(gaze_x^2 + gaze_y^2);
        
        % Check if the new position is outside the boundary of the eyeball
        if distance_from_center > eye_radius
            % Scale the position back so that it lies within the boundary
            scale_factor = eye_radius / distance_from_center;
            gaze_x = gaze_x * scale_factor;
            gaze_y = gaze_y * scale_factor;
        end

        % Clear the previous iris and pupil
        %subplot(1,2,1);
        figure(1)
        cla;
        
        % *************** First Eye (Left Eye) ***************
        % Redraw the first eyeball
        fill(x_eye, y_eye, 'w', 'EdgeColor', 'k');  % White circle with black border
        
        % Calculate new iris position based on the gaze offset for the first eye
        iris_x = iris_radius_x * cos(theta) + gaze_x;
        iris_y = iris_radius_y * sin(theta) + gaze_y;
        
        % Calculate new pupil position based on the gaze offset for the first eye
        % Adjust the size of the pupil based on the scale factor
        pupil_radius_x = default_pupil_radius_x * pupil_scale;
        pupil_radius_y = default_pupil_radius_y * pupil_scale;
        
        pupil_x = pupil_radius_x * cos(theta) + gaze_x;
        pupil_y = pupil_radius_y * sin(theta) + gaze_y;

        % Draw the iris for the first eye
        fill(iris_x, iris_y, iris_color, 'EdgeColor', 'none');  % Brown oval
        
        % Draw the pupil for the first eye with updated size
        fill(pupil_x, pupil_y, pupil_color, 'EdgeColor', 'none');  % Black oval

        % *************** Second Eye (Right Eye) ***************
        % Redraw the second eyeball (shifted by eye_offset)
        fill(x_eye + 2, y_eye, 'w', 'EdgeColor', 'k');  % White circle with black border
        
        % Calculate new iris position for the second eye
        iris_x_2 = iris_radius_x * cos(theta) + gaze_x + 2;  % Shift by 2 for second eye
        iris_y_2 = iris_radius_y * sin(theta) + gaze_y;
        
        % Calculate new pupil position for the second eye
        pupil_x_2 = pupil_radius_x * cos(theta) + gaze_x + 2;
        pupil_y_2 = pupil_radius_y * sin(theta) + gaze_y;

        % Draw the iris for the second eye
        fill(iris_x_2, iris_y_2, iris_color, 'EdgeColor', 'none');  % Brown oval
        
        % Draw the pupil for the second eye with updated size
        fill(pupil_x_2, pupil_y_2, pupil_color, 'EdgeColor', 'none');  % Black oval
        axis tight;

        % ********** Update the Cartesian plot with current gaze **********
        %subplot(1,2,2);
        figure(figure2);
        % Update smooth trajectory in real-time
        if jerk_indices(i)
            % Plot the jerky segment in red
            set(h_gaze, 'XData', gaze_positions(i, 1), 'YData', -1*gaze_positions(i, 2), 'CData', [1 1 1]);
            hold on
        else
            % Plot the smooth segment in blue
            set(h_gaze, 'XData', gaze_positions(i, 1), 'YData', -1*gaze_positions(i, 2), 'CData', [1 1 1]);
            hold on
        end
        
        current_time = i / fps;
        title(sprintf('Time: %.2f s', current_time), 'FontSize', 18);
        
        % If current position is a jerky movement, plot it in red
%         if ismember(i, jerk_indices)
%             h_jerk.XData = [h_jerk.XData, gaze_x];
%             h_jerk.YData = [h_jerk.YData, gaze_y];
%         end

        figure(figure3);

        if i~=1 && i~=length(jerk_indices)
            if jerk_indices(i) %(jerk_indices(i+1)) % %|| %|| (jerk_indices(i-1)) % jerk_indices(i) %
                plot([gaze_positions(i-1, 1) gaze_positions(i, 1)], [-1*gaze_positions(i-1, 2) -1*gaze_positions(i, 2)], 'r-', 'LineWidth', 2);  % Plot the combined x and y trajectory
                hold on;
            else
                plot([gaze_positions(i-1, 1) gaze_positions(i, 1)], [-1*gaze_positions(i-1, 2) -1*gaze_positions(i, 2)], 'b-', 'LineWidth', 2);  % Plot the combined x and y trajectory
                hold on
            end
        end

        set(h_gaze2, 'XData', gaze_positions(i, 1), 'YData', -1*gaze_positions(i, 2), 'CData', [1 1 1]);
        uistack(h_gaze2,'top')
        hold on

        figure(figure4)
        if i~=num_frames && i~=1
            plot(1:i, pupil_sizes(1:i), 'k', 'LineWidth', 2); 
            xlim([0 300])
            ylim([1.1 2])
            xticks([])
            yticks([])
            hold on
        end

        % Capture each frame for the videos
        frame1 = getframe(figure1);
        frame2 = getframe(figure2);
        frame3 = getframe(figure3);
        frame4 = getframe(figure4);
        
        % Write the frames to the video files
        writeVideo(video1, frame1);
        writeVideo(video2, frame2);
        writeVideo(video3, frame3);
        writeVideo(video4, frame4);

        % Pause for animation effect
        pause(0.1);
    end
    
    % Mark the start and end points
    %scatter(gaze_positions(1, 1), gaze_positions(1, 2), 100, 'g', 'filled');  % Start point
    %scatter(gaze_positions(end, 1), gaze_positions(end, 2), 100, 'k', 'filled');  % End point
    
    %legend({'Smooth Gaze Path', 'Jerky Movements', 'Start', 'End'});

close(video1);
close(video2);
close(video3);

end
