function [newMatrix,samplingRateData] = check4interp(matrix, samplingRateData)
   
    samplesPerMilisec = samplingRateData/1000;

    % if downsampling needed, let's just drop every other row.
    if samplesPerMilisec == 1
        newMatrix = matrix; % do nothing

    elseif samplesPerMilisec == 2
        newMatrix = matrix(1:2:end, :);

    elseif samplesPerMilisec == 0.5

        upsampling_factor = 1/samplesPerMilisec; % Change this to downsample

        % Define the interpolation method ('linear', 'nearest', 'spline', 'pchip', etc.)
        % In this case using pchin because it preserves the exact values of the orginal data, 
        % I compared this with 'nearest', 'spline', 'linear' by using actual
        % data group truth and error is mean error is about 0.1 pixels, max
        % error is about 2 pixels and is rare
        % this data will also be smoothed as well in the main script to
        % remove tiny high freq noise
        interpolation_method = 'pchip'; 
        
        % Get the number of rows and columns in the matrix
        [num_rows, num_columns] = size(matrix);
        
        % Define the new number of rows after upsampling or downsampling
        new_num_rows = upsampling_factor * num_rows;
        
        % Create a new row index for the interpolated matrix
        new_row_index = linspace(1, num_rows, new_num_rows);
       
        % Initialize the interpolated matrix
        interpMatrix = zeros(new_num_rows, num_columns);
        
        % Interpolate each row of the matrix
        for col = 1:num_columns
            interpMatrix(:, col) = interp1(1:num_rows, matrix(:, col), new_row_index, interpolation_method);
        end

        % pchip matrix for the x,y interp; & rounding other values (time, pupil) to match precision level of original data.
        newMatrix = [round(interpMatrix(:,1)), interpMatrix(:,2:3), round(interpMatrix(:,4:5))];

    end

    % report new sampling rate
    samplingRateData = samplingRateData/samplesPerMilisec;

end