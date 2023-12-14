rng(0)

% simulated data
noise = randn(1,2000)*5;
ramp0 = repmat(300, [1, 1000])*0.1;
ramp1 = flip(1:300); ramp1 = ramp1*0.1;
ramp2 = zeros(1,500);
ramp3 = 1:200; ramp3 = ramp3*0.1;
ramp = [ramp0, ramp1, ramp2, ramp3]; %, ramp4];
signal = noise + ramp;

figure
subplot(2,2,1)
plot(signal)
hold on
plot(ramp, 'k', 'linewidth', 4)
title('original signal + noise')

y2 = lowpass(signal, 0.001);
stimuluson = 1300; stimulusoff = 1800;
subplot(2,2,2)
plot(y2)
title('STEP 1: lowpass function')

temporalthresh = 70; % ms
stepsize = 1;
onset = 1;

d = nan(1, 2000); % derivative
for i=1:(length(y2)/stepsize)
    try
        temp = mean(diff(y2(onset:onset+temporalthresh)));
        d(i) = temp;
        onset = onset + stepsize;
    catch
        d(i) = d(i);
    end
end

subplot(2,2,3)
% assume max slope occurs between 800 ms before stimuluson
% to stimulus off
[~, idx] = min(d(stimuluson-800:stimulusoff));
idx = idx + stimuluson-800;
% find where slope is 0 after this
temp = d(idx:stimulusoff)>0;
supp = find(temp, 1, 'first') + idx;
%abs_d = abs(d);
%[~, idx2] = min(abs_d(idx:stimulusoff));
%idx2 = idx2 + idx;
scatter(idx, d(idx), 34, 'filled')
hold on
scatter(supp, d(supp), 34, 'filled')
legend('mid-suppression','max suppression','AutoUpdate','off')
hold on
xline(stimuluson)
hold on
xline(stimulusoff)
hold on
plot(1:2000, d)
hold on
yline(0) % no slope
title('STEP 2: steepest slope and max suppresion')

subplot(2,2,4)
plot(y2)
hold on
scatter(idx, y2(idx), 34, 'filled')
hold on
scatter(supp, y2(supp), 34, 'filled')
title('STEP 3: gather original values')

% then you can repeat the steps to find the start of the suppression
% with this array of Y2 -- find best fitting line to those points

% can also repeat these steps to find inflection point

% save, per condition: 
% slope of suppression vector
% start of suprression (ms)
% max suppression (ms)
% point of inflection (ms)



