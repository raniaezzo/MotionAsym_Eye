clear all
clc
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
condition = {'Full_distance_non_radialtangential','Full_distance_radialtangential'};
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';

for ii = 7 :8
    for jj = 2 : 2
        for kk = 1 : 8
            if ii == 3 & jj ==1 & kk ==1 
                continue
            end
                edf2asc = 'F:\RadialBias_pilot1-main\RadialBias_pilot1-main\Analysis\edf2asc';

                main_folder = fullfile('F:\pupildata\Data_DI_wEYE\Data_DI_wEYE', subject{ii}, ...
                    'RawData', condition{jj}, 'Block1');
                csv_filepath = fullfile(main_folder,sprintf('expRes%s_1dMotionAsym_Psychophysics_%s.csv', ...
                    subject{ii}, direction{kk}));
                scr_filepath = fullfile(main_folder,sprintf('scr_file%s_1dMotionAsym_Psychophysics_%s.mat', ...
                    subject{ii}, direction{kk}));
                cd(fullfile(main_folder, 'eyedata')); edf_name = dir(sprintf('*%s*.edf', direction{kk})).name;
                edf_path = fullfile(main_folder,'eyedata',edf_name);
                [path_folder,~,~] = fileparts(edf_path);
                
                MATpath = fullfile(main_folder, 'eyedata','MATs');
                tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
                load(tab_path)
             
                [path, filename, ~] = fileparts(edf_path);
                
                cd(path);
                dat = importdata(sprintf('%s.dat',filename));
                Dat_all_new = []; Dat_stim = [];
                for t = 1:size(tab,1)

                    currTStart = tab(t,2);
                    currTEnd = tab(t,8);
                    currTEnd_est = currTStart + 2200;

                    currGStart = currTStart+1300;
                    currGEnd = currGStart+500;


                    currDat = dat(dat(:,1)>=currTStart & dat(:,1)<=currTEnd_est,:);



                    Dat_all_new = [currDat; Dat_all_new];

                end


                save(sprintf('MATs/%s_Dat_all_new_22.mat',filename),'Dat_all_new');

        end
    end
end