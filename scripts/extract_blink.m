%% this function used to extract blink information used .msg file
function blinkp=extract_blink(filepath)
         [path, filename,~] = fileparts(filepath);
         file=sprintf('%s.msg',filename);
         cd(path);
         msgfid = fopen(file,'r');
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
         blinkp=sprintf('%s\MATs\%s_blink.mat',path,filename);
end