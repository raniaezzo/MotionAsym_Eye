%% create tab file
%%input filepath of eye movement data and the script edf2asc
function [path1,path2,path3]=convert(filename1,filename2)
[path, filename, ~] = fileparts(filename1);
cd(path);
edf2asc=[filename2];
[~,~] = system([edf2asc,' ',filename,'.edf -e -y']); % convert edf to asc
movefile(sprintf('%s.asc',filename), sprintf('%s.msg',filename)); % rename part1 asc to msg (strs)
[~,~] = system([edf2asc,' ',filename,'.edf -s -miss -1.0 -y']);
movefile(sprintf('%s.asc',filename), sprintf('%s.dat',filename)); % rename part2 asc to dat (#s)
msgstr = sprintf('%s.msg',filename);
msgfid = fopen(msgstr,'r');

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
            
            switch char(la(3))
                case 'TRIAL_START'; t = t+1; % count up (add this to the message that is  the first message and occurs reliably in every trial)
                    tab(t,1)  = str2double(char(la(4)));    % trialID
                    tab(t,2)  = str2double(char(la(2)));    % trial start time
                case 'SYNCTIME' % ET: start time: 3150621
                    try
                        tab(t,4)  = str2double(char(la(2)));
                    catch
                        disp('skipping non-trial')
                    end
                case 'STIMULUS_ON' %stimulus on
                    tab(t,5)  = str2double(char(la(2)));
                case 'EVENT_ClearScreen'  %not sure if it is stimulus end
                    tab(t,6)  = str2double(char(la(2)));
                case 'BROKE_FIXATION'
                    tab(t,7)  = str2double(char(la(2)));
                case 'TRIAL_END'
                    tab(t,8)  = str2double(char(la(2)));    
            end
        end
    end
end

% create MAT folder if does not exist
if ~isfolder('MATs')
    mkdir('MATs')
end

save(sprintf('MATs/%s_tab.mat',filename),'tab');
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
         save(sprintf('%s/MATs/%s_blink.mat',path,filename),'blink')
         path3=sprintf('%s/MATs/%s_blink.mat',path,filename);


dat = importdata(sprintf('%s.dat',filename));

tab = tab(tab(:,7)==0,:);
Dat_all = []; Dat_stim = [];
for t = 1:size(tab,1)

    currTStart = tab(t,2);
    currTEnd = tab(t,8);
    
    currGStart = currTStart+1300;
    currGEnd = currGStart+500;
    
    
    currDat = dat(dat(:,1)>=currTStart & dat(:,1)<=currTEnd,:);
    
    currDat_stim = dat(dat(:,1)>=currGStart & dat(:,1)<=currGEnd,:);

    Dat_all = [currDat; Dat_all];
    Dat_stim = [currDat_stim; Dat_stim];
    
%     figure;plot(currDat(:,2),currDat(:,3))
%     figure;plot(currDat_stim(:,2),currDat_stim(:,3))
end

% Delete msg & dat when no more needed (always keep the edf)
%delete(sprintf('%s.msg',filename));
%delete(sprintf('%s.dat',filename));

save(sprintf('MATs/%s_Dat_all.mat',filename),'Dat_all');
save(sprintf('MATs/%s_Dat_stim.mat',filename),'Dat_stim');
path1=sprintf('%s/MATs/%s_tab.mat',path,filename);
path2=sprintf('%s/MATs/%s_Dat_stim.mat',path,filename);
end
