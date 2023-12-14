function convert_new(filename1,filename2)
[path, filename, ~] = fileparts(filename1);
cd(path);
edf2asc=[filename2];
[~,~] = system([edf2asc,' ',filename,'.edf -e -y']); % convert edf to asc
movefile(sprintf('%s.asc',filename), sprintf('%s.msg',filename)); % rename part1 asc to msg (strs)
[~,~] = system([edf2asc,' ',filename,'.edf -s -miss -1.0 -y']);
movefile(sprintf('%s.asc',filename), sprintf('%s.dat',filename)); % rename part2 asc to dat (#s)
msgstr = sprintf('%s.msg',filename);
msgfid = fopen(msgstr,'r');

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
                    tab(t,4)  = str2double(char(la(2)));
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
% save(sprintf('MATs/%s_tab.mat',filename),'tab');
fclose(msgfid);
dat = importdata(sprintf('%s.dat',filename));

tab = tab(tab(:,7)==0,:);
Dat_all_new = []; Dat_stim = [];
for t = 1:size(tab,1)

    currTStart = tab(t,2);
    currTEnd = tab(t,8);
    currTEnd_est = currTStart + 2500;
    
    currGStart = currTStart+1300;
    currGEnd = currGStart+500;
    
    
    currDat = dat(dat(:,1)>=currTStart & dat(:,1)<=currTEnd_est,:);
    
%    currDat_stim = dat(dat(:,1)>=currGStart & dat(:,1)<=currGEnd,:);

    Dat_all_new = [currDat; Dat_all_new];
%    Dat_stim = [currDat_stim; Dat_stim];
    
%     figure;plot(currDat(:,2),currDat(:,3))
%     figure;plot(currDat_stim(:,2),currDat_stim(:,3))
end


save(sprintf('MATs/%s_Dat_all_new.mat',filename),'Dat_all_new');
%save(sprintf('MATs/%s_Dat_stim.mat',filename),'Dat_stim');
end
