clear all; clc; 
subject = {'S01','S02','S03','S04','S05','S06'}; 
condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';
c={'Original','Switched'};
for kk = 3:3
    for ii = 1 : 1
        for jj = 1 : 1
                main_folder = fullfile('F:\RadialBias_pilot1-main\Data_DI_wEYE\Data_DI_wEYE', subject{kk}, ...
                     'RawData', condition{ii}, 'Block1');
                cd(fullfile(main_folder, 'eyedata'));
                edf_name = dir(sprintf('*%s*.edf', direction{jj})).name;
                edf_path = fullfile(main_folder,'eyedata',edf_name);
                msg_filepath=replace(edf_path,'edf','msg');
                samplingRateData=findSamplingRate(msg_filepath);

                MATpath = fullfile(main_folder, 'eyedata','MATs');
        %         cd(fullfile(main_folder, 'eyedata','MATs'));
        %         ms_name = dir(sprintf('%s.mat', direction{j})).name;
                ms_path= fullfile(MATpath,sprintf('%s.mat', direction{jj}) );
                load(ms_path);
                %%
                nTrials = 800;
                nSampleCutOff = 2.5*samplingRateData;
                eventMatrix = nan(nTrials, nSampleCutOff);
                eyetraceMatrix = nan(nTrials, nSampleCutOff,2);

                numBlink=[];

                % figure
                for i = 1 : 800
                    
                    trial_saccades = MS_TEMP(MS_TEMP(:,10)==i, :); 
                    [segs, ~] = size(trial_saccades); 
                    if segs~=0
                        for ti=1:segs
                            saccadeStart = trial_saccades(ti, 8);
                            saccadeEnd = trial_saccades(ti, 9);
                            if saccadeEnd > nSampleCutOff
                                saccadeEnd = nSampleCutOff;
                            end
                            eventMatrix(i, saccadeStart:saccadeEnd) = 2;
                        end
                    end
                end


                EVENTS = eventMatrix(:,1:nSampleCutOff);
                EVENTS(EVENTS==0) = NaN; %blinks
                EVENTS(EVENTS==3) = NaN; %saccades

                EVENTS = EVENTS-1; % 0 to 1
                rate_name = [subject{kk},'_',c{ii},'_',direction{jj},'_rate.mat'];
                rate = nansum(EVENTS,1)./421;
                save(rate_name,'rate')
                figure 
                plot(rate)
                name = [subject{kk},condition{ii},direction{jj}];
                title(name)
        end
    end
end

