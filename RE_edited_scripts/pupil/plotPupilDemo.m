
h = pupilrf(0:4000,10.1,930,0);

figure
%plot(h, 'LineWidth', 5, 'Color', [.5 .5 .5])
shiftal = 500; %1000;
amp = 1.2;
temp = conv(h,amp);
hold on
%plot(temp, 'LineWidth', 5, 'Color', [.3 .3 .3])
hold on
plot([930+shiftal 930+shiftal], [0 amp],':', 'LineWidth', 5, 'Color', [.5 .5 .5])
set(gca, 'FontName', 'Arial', 'FontSize', 12);
ylabel('pupil area (psc)')
xlabel('time (ms)')
xlim([0 4500])
ylim([0, 1.6])
temp = conv(h,amp);
%plot(h, 'LineWidth', 5, 'Color', [.3 .3 .3])
a1 = gca;
a1.XTick = linspace(0, 4500, 4); %[-500 0 500 1000 1500 2000 2500 3000 3500 4000 4500]
%f1 = gcf;
%f1.Position = [217 901 619 354];
plot([shiftal shiftal], [0 amp],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
%hold on
plot([shiftal:length(temp)], [temp(1:end-shiftal+1)], 'LineWidth', 5, 'Color', [.3 .3 .3])
%box on

box off
set(gca, 'XTick', [], 'YTick', [])

% set(gca, 'FontName', 'Arial', 'FontSize', 20);
% set(gca, 'LineWidth', 2);

f1 = gcf;
f1.Position = [82 874 397 285];

% screen_size = get(0, 'ScreenSize');
dpi = get(0, 'ScreenPixelsPerInch');

% set(gcf, 'Position', [1 1 screen_size(3)/2, (screen_size(3)/2)*(2/3)]);
% set(gca, 'FontName', 'Arial', 'FontSize', 64.5);
set(gca, 'LineWidth', 1.5);
% screen_size_inches = (screen_size/2) / dpi;
% %set(gcf, 'PaperPositionMode', 'auto'); 
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 screen_size_inches(3) screen_size_inches(4)]); %[0 0 397 285]);

pathtoSave = '/Users/rje257/Desktop/oculomotor_manuscript_figures/';
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gca, 'FontName', 'Arial', 'FontSize', (285/dpi)*5.6);
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, '-dpdf', fullfile(pathtoSave,'fig3A.pdf'));  % for PDF

%%

figure
% boxcarEvent = zeros(length(5000), 1);
% boxcarEvent(1:500) = .01; % 1%
% temp = conv(h,boxcarEvent);
% hold on
% plot([zeros(1,500), temp(1:end-500)], 'LineWidth', 5, 'Color', [.3 .3 .3])
x_patch = [0, 500, 500, 0]; % x-coordinates
y_patch = [0, 0, max(temp), max(temp)]; % y-coordinates
fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
hold on
x_patch = [-1300, 0, 0, -1300]; % x-coordinates
y_patch = [0, 0, max(temp), max(temp)]; % y-coordinates
fill(x_patch, y_patch, [255, 255, 0]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
x_patch = [500, 4500, 4500, 500]; % x-coordinates
y_patch = [0, 0, max(temp), max(temp)]; % y-coordinates
fill(x_patch, y_patch, [255, 0, 0]/255, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
hold on
plot([1500 1500], [0 1],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
hold on
plot([250 250], [0 1],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
hold on
plot([-300 -300], [0 1],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
hold on
plot([0 0], [0 .5],'-', 'LineWidth', 5, 'Color', [0.3 .3 .3])
hold on
plot([500 500], [0 .5],'-', 'LineWidth', 5, 'Color', [0.3 .3 .3])
hold on
plot([-15 515], [.5 .5],'-', 'LineWidth', 5, 'Color', [0.3 .3 .3])

ylabel('event impulses')
xlabel('time (ms)')
xlim([-1300 2500])
xlim([-500 2500])
a1 = gca;
a1.XTick = -1000:500:2500; 
f1 = gcf;
% f1.Position = [217 901 619 354];
f1.Position = [73 983 795 285];

ylim([0 1.2])
yticks([]);
% set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
box off
set(gca, 'YTick', [])

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gca, 'FontName', 'Arial', 'FontSize', (285/dpi)*5.6);
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, '-dpdf', fullfile(pathtoSave,'fig3B.pdf'));  % for PDF

%% load

load('/Volumes/Vision/UsersShare/Rania/MS_Project/S01trial35_examplePupilFit.mat')
load('/Volumes/Vision/UsersShare/Rania/MS_Project/Data_DI_wEYE/S01/ProcessedData/Summary/pupil/S01_allpupilData.mat')
figure
sj = pret_estimate_sj(sj,model,wnum,options);
xline(800, ':', 'Color', [1 0 0 0.05], 'LineWidth', 3)
hold on
x_patch = [0, 500, 500, 0]; % x-coordinates
y_patch = [0, 0, 200, 200]; % y-coordinates
fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
ylim([0 200])
hold on
legend({'data', '', '', '', '', '', 'model'});
legend box off
f1 = gcf;
f1.Position = [53 571 795-(795*.21) 285];

% correction of psc
yticks = get(gca, 'YTick');
% Set the Y-tick labels to be divided by 10
yticks = linspace(0,20,5) / (mean(alltrialSignalSummary(:,2))/max(alltrialSignalSummary(:,4)));
set(gca, 'YTickLabel', yticks * (mean(alltrialSignalSummary(:,2))/max(alltrialSignalSummary(:,4)))); % / 10);

% set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
box off
xlabel('time (ms)')
ylabel('pupil area (psc)')

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gca, 'FontName', 'Arial', 'FontSize', (285/dpi)*5.6);
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, '-dpdf', fullfile(pathtoSave,'fig3C_new.pdf'));  % for PDF