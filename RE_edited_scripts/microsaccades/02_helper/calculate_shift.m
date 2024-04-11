%% This script is not used - but was to explore different ways to capture the lag.

% Assuming data1 and data2 are your time series

data1 = nanmean(rate.cardinal,1)*60;
data2 = nanmean(rate.oblique,1)*60;

% data1 = data1(300:1000);
% data2 = data2(300:1000);

data1 = data1(550:1300);
data2 = data2(550:1300);

% goal: see if a stretch or a shift better explains the data transofrmation
% from data1 to data2

figure
plot(data1)
[val,idx] = min(data1(300:800));
hold on
xline(idx+300)
hold on
xline(300, '--')
hold on
xline(800, '--')

% event 1 : -300 when the fixation is enforced
% event 2 : 0 when stimulus is on

%% find the correlation between arrays when shifting in time

% Compute cross-correlation
[corr_values, lag] = xcorr(data1, data2);

% Find lag at which cross-correlation is maximized
[~, max_index] = max(abs(corr_values));
optimal_shift = lag(max_index);

% Shift data2 by the optimal shift
data2_shifted = circshift(data2, [0, optimal_shift]);

% Now data1 and data2_shifted are aligned
figure;
subplot(1,2,1)
plot(data1, 'g', 'Linewidth', 2)
hold on
plot(data2, 'b', 'Linewidth', 2)
hold on
plot(data2_shifted, 'b--', 'Linewidth', 2)

%% another way to compute shift

% Define a range of shifts to consider
max_shift = 100; %min(length(data1), length(data2)) - 1; % Maximum shift allowed
shifts = -max_shift:max_shift;

% Compute the differences for each shift
differences = zeros(size(shifts));
for i = 1:length(shifts)
    shift = shifts(i);
    if shift >= 0
        differences(i) = sum(abs(data1(1:end-shift) - data2(shift+1:end)));
    else
        differences(i) = sum(abs(data1(-shift+1:end) - data2(1:end+shift)));
    end
end

% Find the optimal shift that minimizes the difference
[~, min_index] = min(differences);
optimal_shift = shifts(min_index);

% Shift data2 by the optimal shift
if optimal_shift >= 0
    data2_shifted = [zeros(1, optimal_shift), data2(1:end-optimal_shift)];
else
    data2_shifted = [data2(-optimal_shift+1:end), zeros(1, -optimal_shift)];
end

% Now data1 and data2_shifted are aligned
figure;
subplot(1,2,1)
plot(data1, 'g', 'Linewidth', 2)
hold on
plot(data2, 'b', 'Linewidth', 2)
hold on
plot(data2_shifted, 'b--', 'Linewidth', 2)

%% scaling factor that maximizes the similarity between the two signals after stretching or compressing one of them

% Compute cross-correlation
[corr_values, scale_factors] = xcorr(data1, data2);

% Find scale factor at which cross-correlation is maximized
[~, max_index] = max(abs(corr_values));
optimal_scaling_factor = scale_factors(max_index);

% Display result
fprintf('Optimal Time Scaling Factor: %.4f\n', optimal_scaling_factor);


%% Additive Shift

% Define the cost function
cost_function = @(params) norm(data1 - circshift(data2, round(params(1))));

% Initial guess for the parameter (temporal shift)
initial_guess = -25;  % Initial guess for the temporal shift

% Define lower and upper bounds for the parameter
lb = -100; %length(data2);  % Lower bound
ub = 100; %length(data2);   % Upper bound

% Options for optimization (optional)
options = optimset('Display','iter');  % Display optimization process

% Run optimization
optimal_temporal_shift = fmincon(cost_function, initial_guess, [], [], [], [], lb, ub, [], options);

% Apply the estimated temporal shift to data2
data2_fit = circshift(data2, round(optimal_temporal_shift));

% Compute the Euclidean distance between data1 and data2_fit
euclidean_distance = norm(data1 - data2_fit);

% Display results
fprintf('Optimal Temporal Shift: %.4f\n', optimal_temporal_shift);
fprintf('Euclidean Distance: %.4f\n', euclidean_distance);


figure;
subplot(1,2,1)
% plot(lag, corr_values);
% xlabel('Lag');
% ylabel('Cross-correlation');
% title('Cross-correlation Function');
% grid on;
% subplot(1,2,2)
plot(data1, 'g', 'Linewidth', 2)
hold on
plot(data2, 'b', 'Linewidth', 2)
hold on
plot(data2_fit, 'b--', 'Linewidth', 2)




%% Multiplicative factor & shift

% Define the cost function
cost_function = @(params) norm(data1 - params(1) * circshift(data2, round(params(2))));

% Initial guesses for parameters (multiplicative factor and horizontal shift)
initial_params = [1, 0];  % [factor, shift]

% Run optimization
optimal_params = fminsearch(cost_function, initial_params);

% Extract estimated parameters
multiplicative_factor = optimal_params(1);
horizontal_shift = round(optimal_params(2)); % Round to nearest integer

% Apply estimated parameters to data2
data2_fit = multiplicative_factor * circshift(data2, horizontal_shift);

% Compute the Euclidean distance between data1 and data2_fit
euclidean_distance = norm(data1 - data2_fit);

% Display results
fprintf('Estimated Multiplicative Factor: %.4f\n', multiplicative_factor);
fprintf('Estimated Horizontal Shift: %.4f\n', horizontal_shift);
fprintf('Euclidean Distance: %.4f\n', euclidean_distance);

figure;
subplot(1,2,1)
% plot(lag, corr_values);
% xlabel('Lag');
% ylabel('Cross-correlation');
% title('Cross-correlation Function');
% grid on;
% subplot(1,2,2)
plot(data1, 'g', 'Linewidth', 2)
hold on
plot(data2, 'b', 'Linewidth', 2)
hold on
plot(data2_fit, 'b--', 'Linewidth', 2)

%%
% 
% % Normalize data
% data1 = (data1 - mean(data1)) / std(data1);
% data2 = (data2 - mean(data2)) / std(data2);
% 
% % Compute cross-correlation
% [corr_values, lag] = xcorr(data1, data2);
% 
% % Find peak of cross-correlation
% [max_corr, max_index] = max(abs(corr_values));
% 
% % Estimate temporal shift
% temporal_shift = lag(max_index);
% 
% % Estimate amplitude shift (optional)
% amplitude_shift = max_corr;
% 
% % Display results
% fprintf('Temporal Shift: %d\n', temporal_shift);
% fprintf('Amplitude Shift: %.4f\n', amplitude_shift);
% 
% % Plot cross-correlation function
% figure;
% plot(lag, corr_values);
% xlabel('Lag');
% ylabel('Cross-correlation');
% title('Cross-correlation Function');
% grid on;

%%
% Assuming data1 and data2 are your time series

% Normalize data
data1 = (data1 - mean(data1)) / std(data1);
data2 = (data2 - mean(data2)) / std(data2);

% Compute cross-correlation
[corr_values, lag] = xcorr(data1, data2);

% Find peak of cross-correlation
[~, max_index] = max(abs(corr_values));
temporal_shift = lag(max_index);

% Apply temporal shift to data2
data2_shifted = circshift(data2, temporal_shift);

% Trim data1 and data2_shifted to match lengths
min_length = min(length(data1), length(data2_shifted));
data1_trimmed = data1(1:min_length);
data2_shifted_trimmed = data2_shifted(1:min_length);

% Perform linear regression
X = [ones(length(data1_trimmed), 1), data1_trimmed(:)]; % Transpose data1_trimmed if necessary
coefficients = X\data2_shifted_trimmed(:);
amplitude_shift = coefficients(2);  % Slope coefficient represents amplitude shift

% Display results
fprintf('Temporal Shift: %d\n', temporal_shift);
fprintf('Amplitude Shift: %.4f\n', amplitude_shift);

%%
% Plot cross-correlation function (optional)
figure;
subplot(1,2,1)
plot(lag, corr_values);
xlabel('Lag');
ylabel('Cross-correlation');
title('Cross-correlation Function');
grid on;
subplot(1,2,2)
plot(data1, 'g', 'Linewidth', 2)
hold on
plot(data2, 'b', 'Linewidth', 2)
hold on
plot(data2_shifted, 'b--', 'Linewidth', 2)