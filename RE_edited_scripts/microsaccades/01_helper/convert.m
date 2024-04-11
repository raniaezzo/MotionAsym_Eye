%% create tab file
%%input filepath of eye movement data and the script edf2asc
function [path1,path2,path3,path4,samplingRateData]=convert(filename1,filename2,savetabfolder)
    [path, filename, ~] = fileparts(filename1);
    cd(path);
    edf2asc=[filename2];
    [~,~] = system([edf2asc,' ',filename,'.edf -e -y']); % convert edf to asc
    movefile(sprintf('%s.asc',filename), sprintf('%s.msg',filename)); % rename part1 asc to msg (strs)
    [~,~] = system([edf2asc,' ',filename,'.edf -s -miss -1.0 -y']);
    movefile(sprintf('%s.asc',filename), sprintf('%s.dat',filename)); % rename part2 asc to dat (#s)
    msgstr = sprintf('%s.msg',filename);
    msgfid = fopen(msgstr,'r');

    samplingRateData=findSamplingRate(msgstr); % get sampling rate
    samplesPerMilisec = samplingRateData/1000;
    
    % path_blink=sprintf('%s\%s',path,msgstr);
    % path3=extract_blink(path_blink);
    
    tab = [];
    t = 0; % Rania: 0
    
    sameData = 1;
    while sameData
        line = fgetl(msgfid);
        if ~ischar(line)                            % end of file
            sameData = 0;
            break;
        end
        if ~isempty(line)                           % skip empty lines
            la = strread(line,'%s');                % matrix of strings in line
            if length(la) >= 3
                
                %char(la(3))
                
                charArray = char(la(3));
                charTimeStamp = char(la(2));

                %switch char(la(3))
                    %case 'TRIAL_START'; t = t+1; % count up (add this to the message that is  the first message and occurs reliably in every trial)
                    if strcmp(charArray, 'TRIAL_START')
                        t = t+1;
                        tab(t,1)  = str2double(char(la(4)));    % trialID
                        tab(t,2)  = str2double(charTimeStamp);    % trial start time
                    %case 'SYNCTIME' % ET: start time: 3150621
                    elseif strcmp(charArray, 'SYNCTIME')
                        try
                            tab(t,4)  = str2double(charTimeStamp);
                        catch
                            disp('skipping non-trial')
                        end
                    %case 'STIMULUS_ON'
                    elseif strcmp(charArray, 'STIMULUS_ON')
                        tab(t,5)  = str2double(charTimeStamp);
                    %case 'EVENT_ClearScreen' 
                    elseif strcmp(charArray, 'EVENT_ClearScreen')
                        tab(t,6)  = str2double(charTimeStamp);
                    %case 'BROKE_FIXATION'
                    elseif strcmp(charArray, 'BROKE_FIXATION')
                        tab(t,7)  = str2double(charTimeStamp);
                    %case 'TRIAL_END'
                    elseif strcmp(charArray, 'TRIAL_END')
                        tab(t,8)  = str2double(charTimeStamp);    
                    end
            end
        end
    end

    % FIX (WORK AROUND): ensure that when Eyelink channel receive error occurs, this does
    % not mess up the trials. I can do tab(t,8) - tab(t, 2) to double check
    % the duration per trial is reasonable and NOT negative or 0.

    fixInds = tab(tab(:,8) - tab(:,2) < 0);

    for fi=1:length(fixInds)
        idx = fixInds(fi);
        if tab(idx,8) == 0
            tab(idx,8) = tab(idx+1,2); % start of the next trial
            warning(sprintf('Eyelink channel clogged for trial %i. Changing TRIAL_END=0 to next TRIAL_START.', fi))
        elseif tab(idx,8)<tab(idx,2)
            tab(idx,8) = tab(idx+1,2); % start of the next trial
            warning(sprintf('Eyelink channel clogged for trial %i. Setting TRIAL_END<0 to next TRIAL_START.', fi))
        else
            error('Unencountered issue with the ASC / MSG conversion. Check convert.m for details..')
        end
    end
    
    
%     save(sprintf('MATs/%s_tab.mat',filename),'tab');
    path1=fullfile(savetabfolder, sprintf('%s_tab.mat',filename));
    save(path1,'tab')
    fclose(msgfid);
    
    % path_blink=sprintf('%s\%s',path,msgstr);
    % path3=extract_blink(path_blink);
    msgfid = fopen(msgstr,'r');
             blink = [];
             t = 0; % Rania: 0
             sameData = 1;
             while sameData
                    line = fgetl(msgfid);
                    if ~ischar(line)                            % end of file
                        sameData = 0;
                        break;
                    end
                    if ~isempty(line)                           % skip empty lines
                        la = strread(line,'%s');                % matrix of strings in line
                        if length(la) >= 3
                            switch char(la(1))
                                case 'EBLINK'; t = t+1; % count up (add this to the message that is  the first message and occurs reliably in every trial)
                                    blink(t,1)  = str2double(char(la(3)));   
                                    blink(t,2)  = str2double(char(la(4)));
                                    blink(t,3)  = str2double(char(la(5)));
                         
                            end
                        end
                    end
             end
             path2 = fullfile(savetabfolder, sprintf('%s_blink.mat', filename));
             save(path2, 'blink')
    
    
    dat = importdata(sprintf('%s.dat',filename));

    % in terms of samples
    stimulusON = 1300*samplesPerMilisec; % from trial start
    stimuluslength = 500*samplesPerMilisec;

    tab = tab(tab(:,7)==0,:);
    Dat_all = []; Dat_stim = [];
    for t = 1:size(tab,1)
    
        currTStart = tab(t,2);
        currTEnd = tab(t,8);
        
        currGStart = currTStart+stimulusON;
        currGEnd = currGStart+stimuluslength;
        
        
        currDat = dat(dat(:,1)>=currTStart & dat(:,1)<=currTEnd,:);
        
        currDat_stim = dat(dat(:,1)>=currGStart & dat(:,1)<=currGEnd,:);
    
        Dat_all = [Dat_all; currDat];
        Dat_stim = [Dat_stim; currDat_stim];

%         Dat_all = [currDat; Dat_all];
%         Dat_stim = [currDat_stim; Dat_stim];
        
    %     figure;plot(currDat(:,2),currDat(:,3))
    %     figure;plot(currDat_stim(:,2),currDat_stim(:,3))
    end
    
    
    path3=fullfile(savetabfolder, sprintf('%s_Dat_stim.mat',filename));
    save(path3,'Dat_stim');

    path4=fullfile(savetabfolder, sprintf('%s_Dat_all.mat',filename));
    save(path4,'Dat_all');

end
