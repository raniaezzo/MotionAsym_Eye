load('/Users/rje257/Desktop/transfer_gate/S03/S03_allmsData.mat')

% Create example binary matrix (0's and 1's)
matrix = allMSData;  % Example 10x10 matrix with random 0's and 1's
matrix(isnan(matrix)) = 0;
[nrows, ncols] = size(matrix);
matrix = logical(matrix);
matrix = matrix(1:5500,500:2250);

% Create example binary matrix (0's and 1's)
%matrix = randi([0, 1], 10, 10);  % Example 10x10 matrix with random 0's and 1's

% Set up figure for displaying the frames
figure;
ylim([0 nrows])
xlim([0 ncols])

% Number of frames for the GIF
num_frames = nrows;

% Initialize matrix for falling simulation
falling_matrix = matrix;

% GIF settings
filename = 'falling_ones.gif';

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
    colormap(gray);  % Set colormap to grayscale for better visibility
    %axis equal;
    axis off;
    title(['Frame ' num2str(t)]);
    drawnow;
    
    % Capture the frame for the GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if t == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end

% Display final "pile" in histogram
pile_heights = sum(falling_matrix == 1, 1);
figure;
bar(pile_heights);
xlabel('Column');
ylabel('Pile Height');
title('Final Pile Height Histogram');

%%

% Load the GIF
[gif, cmap] = imread(filename, 'Frames', 'all');

% Create VideoWriter object for MP4 file
v = VideoWriter('new_falling_ones.mp4', 'MPEG-4');
v.FrameRate = 100 * 10;  % Original frame rate * 100 for 100x speedup
open(v);

% Write each frame to the video
num_frames = size(gif, 4);
for k = 1:num_frames
    rgb_frame = ind2rgb(gif(:, :, 1, k), cmap);  % Convert indexed to RGB
    writeVideo(v, rgb_frame);
end

% Close the video file
close(v);