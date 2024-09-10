





[Ycalc, X] = pret_calc(testing(1)); 

figure
plot(Ycalc(1300:end))
%%

pret_plot_model(testing(1))

%%

close all;
% create an ave_struct
% easydata = pupilFit_easy(1);

rawEasy = pupilData(easyTrials,:);
rawHard = pupilData(hardTrials,:);

varExp_easy = cell2mat({pupilFit_easy.('R2')});
varExp_hard = cell2mat({pupilFit_hard.('R2')});

[easyVal, easyIdx] = max(varExp_easy);
[hardVal, hardIdx] = max(varExp_hard);

easyIdx = 5;
hardIdx = 3;

easydata = pupilFit_easy(easyIdx);
harddata = pupilFit_hard(hardIdx);

easyTrialraw = rawEasy(easyIdx,:);
hardTrialraw = rawHard(hardIdx,:);

% pupilFit_easy_filt = pupilFit_easy(varExp_easy>0.5);
% easydata.eventtimes = median(cell2mat({pupilFit_easy_filt.('eventtimes')}'), 1);
% easydata.ampvals = median(cell2mat({pupilFit_easy_filt.('ampvals')}'), 1);
% easydata.boxampvals = median(cell2mat({pupilFit_easy_filt.('boxampvals')}'), 1);
% easydata.latvals = median(cell2mat({pupilFit_easy_filt.('latvals')}'), 1);
% easydata.R2 = median(varExp_easy);
% easydata.cost = []; easydata.BICrel = [];


pret_plot_model(easydata)
hold on
plot(easyTrialraw(1300:end), 'Color', 'g', 'LineWidth', 2)
%plot(nanmedian(pupilData(easyTrials, 1300:end),1), 'Color', 'g', 'LineWidth', 2)
% figure
% pret_plot_model(harddata)
% hold on
% plot(nanmedian(pupilData(hardTrials, 1300:end),1), 'Color', 'r', 'LineWidth', 2)

figure
pret_plot_model(harddata)
hold on
plot(hardTrialraw(1300:end), 'Color', 'r', 'LineWidth', 2)