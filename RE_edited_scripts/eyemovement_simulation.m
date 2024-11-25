% Example gaze positions (moving from center to top-right)
gaze_positions = [0, 0; 0.1, 0.1; 0.2, 0.2; 0.3, 0.3; 0.4, 0.4];

% Example pupil sizes (alternating larger and smaller)
pupil_sizes = [1; 1.2; 0.8; 1.3; 0.9];

% Call the function
rotating_eye(gaze_positions, pupil_sizes);

function rotating_eye(gaze_positions, pupil_sizes)
    % Function to display a cartoon eye that rotates based on gaze positions
    % and adjusts pupil size
    % INPUT:
    % gaze_positions - Nx2 matrix where each row is [x, y] representing the gaze position
    % relative to the center of the screen (normalized between -1 and 1)
    % pupil_sizes - Nx1 vector where each element represents the size (scale) of the pupil

    % Create the figure window
    figure;
    axis equal;
    axis([-1.5 1.5 -1.5 1.5]);  % Set axis limits
    hold on;
    axis off;  % Hide the axis

    % Eyeball (white circle)
    theta = linspace(0, 2*pi, 100);
    x_eye = cos(theta);
    y_eye = sin(theta);
    fill(x_eye, y_eye, 'w', 'EdgeColor', 'k');  % White circle with black border

    % Iris (brown oval)
    iris_radius_x = 0.5; %0.4;
    iris_radius_y = 0.5;
    iris_color = [0.6 0.3 0];  % Brownish color

    % Default Pupil (black oval) size
    default_pupil_radius_x = 0.15; %0.1;
    default_pupil_radius_y = 0.15;
    pupil_color = 'k';  % Black color

    % Loop through each gaze position and update the iris and pupil positions
    for i = 1:size(gaze_positions, 1)
        % Get the current gaze position
        gaze_x = gaze_positions(i, 1);  % Gaze position on x-axis (normalized)
        gaze_y = gaze_positions(i, 2);  % Gaze position on y-axis (normalized)
        
        % Get the current pupil size
        pupil_scale = pupil_sizes(i);  % Scale factor for pupil size
        
        % Scale the gaze position to control the iris and pupil movement
        iris_offset_x = 0.5 * gaze_x;  
        iris_offset_y = 0.5 * gaze_y;
        
        % Clear the previous iris and pupil
        cla;
        
        % Redraw the eyeball
        fill(x_eye, y_eye, 'w', 'EdgeColor', 'k');  % White circle with black border
        
        % Calculate new iris position based on the gaze offset
        iris_x = iris_radius_x * cos(theta) + iris_offset_x;
        iris_y = iris_radius_y * sin(theta) + iris_offset_y;
        
        % Calculate new pupil position based on the gaze offset
        % Adjust the size of the pupil based on the scale factor
        pupil_radius_x = default_pupil_radius_x * pupil_scale;
        pupil_radius_y = default_pupil_radius_y * pupil_scale;
        
        pupil_x = pupil_radius_x * cos(theta) + iris_offset_x;
        pupil_y = pupil_radius_y * sin(theta) + iris_offset_y;

        % Draw the iris
        fill(iris_x, iris_y, iris_color, 'EdgeColor', 'none');  % Brown oval
        
        % Draw the pupil with updated size
        fill(pupil_x, pupil_y, pupil_color, 'EdgeColor', 'none');  % Black oval

        % Pause for animation effect
        pause(0.1);
    end
end