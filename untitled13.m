%% Burlingham/Heeger model:

% Define the dimensions
num_rows = 800;
num_cols = 2;
spacing = 4500;

% Initialize the matrix
result_matrix = zeros(num_rows, num_cols);

% Generate the matrix
for i = 1:num_rows
    % Calculate the starting index for the current row
    start_index = (i - 1) * spacing + 1;
    
    % Assign values to the current row
    result_matrix(i, 1) = start_index;
    result_matrix(i, 2) = start_index + spacing - 1;
end

% Display the result
disp(result_matrix);

%%
VUp = load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Full_distance_non_radialtangential/eyedata/components/S01VU_pupilMatrix.mat');
ULp = load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Full_distance_non_radialtangential/eyedata/components/S01UL_pupilMatrix.mat');

VUeye = load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Full_distance_non_radialtangential/eyedata/components/S01VU_eyetraceMatrix.mat');
ULeye = load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Full_distance_non_radialtangential/eyedata/components/S01UL_eyetraceMatrix.mat');
VUeye_x = VUeye.eyetraceMatrix(:,:,1);
VUeye_y = VUeye.eyetraceMatrix(:,:,2);
ULeye_x = ULeye.eyetraceMatrix(:,:,1);
ULeye_y = ULeye.eyetraceMatrix(:,:,2);

VUtab = load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Full_distance_non_radialtangential/eyedata/MATs/VU080715_tab.mat');
ULtab = load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Full_distance_non_radialtangential/eyedata/MATs/UL080614_tab.mat');

VUtrials = [VUtab.tab(:,2)+1300, VUtab.tab(:,2)+2500];
ULtrials = [ULtab.tab(:,2)+1300, ULtab.tab(:,2)+2500];

in2.xPos = {reshape(VUeye_x, 1.', [])', reshape(ULeye_x, 1.', [])'};
in2.yPos = {reshape(VUeye_y, 1.', [])', reshape(ULeye_y, 1.', [])'};
in2.trialTypes = {ones(1,800), ones(1,800)*2};
in2.startInds = {result_matrix, result_matrix}; %{VUtrials, ULtrials};
in2.pupilArea = {reshape(VUp.pupilMatrix.', 1, [])', reshape(ULp.pupilMatrix.', 1, [])'};
in2.sampleRate = {1000,1000};

%%

fitModel(in2) % fits do not really make sense
%%
figure
plot(VUp.pupilMatrix(1,1:4000), 'r')
% tmp = reshape(VUp.pupilMatrix.', 1, []);
% size(tmp)
% figure
% plot(tmp(1:4000))
hold on
plot(in2.pupilArea{1}(1:4000), 'k--')