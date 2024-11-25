angles = linspace(0, 360, 100);  % Rotate from 0 to 360 degrees
rotating_eye(angles);

function rotating_eye(angles)
    % Function to display a rotating cartoon eye based on angular values
    % INPUT:
    % angles - vector of angular values (in degrees) for the pupil/iris rotation
    
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
    iris_radius_x = 0.4;
    iris_radius_y = 0.5;
    iris_color = [0.6 0.3 0];  % Brownish color

    % Pupil (black oval)
    pupil_radius_x = 0.1;
    pupil_radius_y = 0.15;
    pupil_color = 'k';  % Black color

    % Loop through each angle value and rotate the iris and pupil
    for angle = angles
        % Convert angle to radians
        angle_rad = deg2rad(angle);
        
        % Clear the previous iris and pupil
        cla;
        
        % Redraw the eyeball
        fill(x_eye, y_eye, 'w', 'EdgeColor', 'k');  % White circle with black border
        
        % Calculate new iris and pupil positions based on the angle
        iris_x = iris_radius_x * cos(theta) + 0.5 * cos(angle_rad);
        iris_y = iris_radius_y * sin(theta) + 0.5 * sin(angle_rad);
        
        pupil_x = pupil_radius_x * cos(theta) + 0.5 * cos(angle_rad);
        pupil_y = pupil_radius_y * sin(theta) + 0.5 * sin(angle_rad);

        % Draw the iris
        fill(iris_x, iris_y, iris_color, 'EdgeColor', 'none');  % Brown oval
        
        % Draw the pupil
        fill(pupil_x, pupil_y, pupil_color, 'EdgeColor', 'none');  % Black oval

        % Pause for animation effect
        pause(0.1);
    end
end
