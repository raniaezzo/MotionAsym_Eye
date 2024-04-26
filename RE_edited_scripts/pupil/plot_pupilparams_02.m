clc; clear all; close all;

% replace path to wherever json file lives
cd('/Users/rje257/Documents/GitHub/MotionAsym_Eye/RE_edited_scripts/')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 

varExpThresh = 0.5;
difficulty = 'direction'; %'tilt'; %'tilt'; % direction
event_label = {'pre_stim', 'sti on', 'response'};

%% Plot the amplitude and latency for each event

for si=1:1 %length(subjects)
    subj=subjects{si};

    savePath = sprintf('%s/%s/ProcessedData/Summary/pupil', datadir, subj);
    load(fullfile(savePath,sprintf('%s_allpupilFits.mat', subj)), 'output');
    load(fullfile(savePath,sprintf('%s_allpupilData.mat', subj)));
    
    % this is the order of the values
    locations = 0:45:315; 
    directions = 0:90:270; 
    locationdegreesArray = [315,135,225,45,270,90,180,0];
    directiondegreesArray = [45,225,135,315,90,270,0,180];
    
    locationIndices = find(ismember(locationdegreesArray, locations));
    carddirectionIndices = find(ismember(directiondegreesArray, directions));
    
    if strcmp(difficulty, 'tilt')
            easyTrials = alltab(:,11)>=4 & alltab(:,14)>=1; % for tilt (correct only)
            hardTrials = alltab(:,11)<8 & alltab(:,14)>=1; % for tilt (correct only)
    elseif strcmp(difficulty, 'direction')
            easyTrials = ismember(alltab(:,10), carddirectionIndices) & alltab(:,14)>=1; %[alltab(2:end,14)>=1; 0]; % for card (correct only)
            hardTrials = ~ismember(alltab(:,10), carddirectionIndices) & alltab(:,14)>=1; %[alltab(2:end,14)>=1; 0]; % for card (correct only)
    elseif strcmp(difficulty, 'outcome')
            easyTrials = [alltab(2:end,14)>=1; 0]; %alltab(:,14)>=1; % for correct
            hardTrials = [alltab(2:end,14)>=0; 0]; %alltab(:,14)>=0;
    end
    
    varExp = {output.('R2')};
    varExpMat = vertcat(varExp{:});
    varExpFilter = varExpMat>varExpThresh; %5;
    
    easyYint = alltrialSignalSummary(easyTrials & varExpFilter, 1);
    hardYint = alltrialSignalSummary(hardTrials & varExpFilter, 1);
    
    jitterEasy = (rand(length(easyYint), 1)-0.5)*.2; jitterHard = (rand(length(hardYint), 1)-0.5)*.2;
    
    figure
    scatter(ones(length(easyYint), 1)+jitterEasy, easyYint(:,1), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
    hold on
    scatter((ones(length(hardYint), 1)*2)+jitterHard, hardYint(:,1), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
    xlim([0 3])
    
    %
    plotOn=1;
    
    if plotOn
       
        
        events = {output.('eventtimes')};
        eventMat = vertcat(events{:});
        tmaxvals = {output.('tmaxval')};
        tmaxMat = vertcat(tmaxvals{:});
    
        amps = {output.('ampvals')};
        ampMat = vertcat(amps{:});
        
        easyAmp = ampMat(easyTrials & varExpFilter, :);
        hardAmp = ampMat(hardTrials & varExpFilter, :);
        
        latencies = {output.('latvals')};
        latMat = vertcat(latencies{:});
        
        easyLat = latMat(easyTrials & varExpFilter, :);
        hardLat = latMat(hardTrials & varExpFilter, :);
        
        boxamp = {output.('boxampvals')};
        boxMat = vertcat(boxamp{:});
        
        easyBox = boxMat(easyTrials & varExpFilter, :);
        hardBox = boxMat(hardTrials & varExpFilter, :);
    
        tmaxEasy = tmaxMat(easyTrials & varExpFilter, :);
        tmaxHard = tmaxMat(hardTrials & varExpFilter, :);
    
        eventsEasy = eventMat(easyTrials & varExpFilter, :);
        eventsHard = eventMat(hardTrials & varExpFilter, :);
        
        %close all;
        figure
        ploti= 1;
        for ee=1:length(event_label)
            subplot(2, length(event_label)+1,ploti)
            %scatter([ones(length(easyAmp), 1), ones(length(hardAmp),1)*2], [easyAmp(:,1), hardAmp(:,1)])
            
            scatter(ones(length(easyAmp), 1)+jitterEasy, easyAmp(:,ee), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
            hold on
            scatter(1, median(easyAmp(:,ee),1), 150, 'ok', 'filled')
            hold on
            scatter((ones(length(hardAmp), 1)*2)+jitterHard, hardAmp(:,ee), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
            hold on
            scatter(2, median(hardAmp(:,ee),1), 150, 'ok', 'filled')
            xlim([0 3])
            ylim([0 150])
            title(sprintf('Amplitude event #%i',ee))
            ylabel('PSC')
            hold off
        
            ploti=ploti+1;
        end
        ee=1;
        subplot(2, length(event_label)+1,ploti)
        scatter(ones(length(easyBox), 1)+jitterEasy, easyBox(:,ee), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
        hold on
        scatter(1, mean(easyBox(:,ee),1), 150, 'ok', 'filled')
        hold on
        scatter((ones(length(hardBox), 1)*2)+jitterHard, hardBox(:,ee), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
        hold on
        scatter(2, mean(hardBox(:,ee),1), 150, 'ok', 'filled')
        ylim([0 100])
        xlim([0 3])
        title(sprintf('Box Amp #%i',ee))
        ylabel('PSC')
        hold off
        ploti=ploti+1;
        for ee=1:length(event_label)
            subplot(2, length(event_label)+1,ploti)
        
            % latency
            scatter(ones(length(easyLat), 1)+jitterEasy, tmaxEasy+easyLat(:,ee)+eventsEasy(:,ee), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
            hold on
            scatter(1, median(tmaxHard)+median(easyLat(:,ee),1)+median(eventsEasy(:,ee)), 150, 'ok', 'filled')
            hold on
            scatter((ones(length(hardLat), 1)*2)+jitterHard, tmaxHard+hardLat(:,ee)+eventsHard(:,ee), 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
            hold on
            scatter(2, median(tmaxHard)+median(hardLat(:,ee),1)+median(eventsHard(:,ee)), 150, 'ok', 'filled')
            xlim([0 3])
            if ee == 1
                ylim([0 2000+max(eventsHard(:,ee))]) % 0 is model start - just use hard
            elseif ee ==2
                ylim([min(tmaxHard) 2000+max(eventsHard(:,ee))]) 
            end
            ploti=ploti+1;
            title(sprintf('Peak event #%i',ee))
            ylabel('ms')
            hold off
        end
        
    
    
    f1 = gcf;
    f1.Position = [24 219 1433 1117];
    sgtitle(difficulty)
    
    end
end