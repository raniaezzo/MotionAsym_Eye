%this script is for the pupil size data process.
clc; clear all;
subject = {'S01','S02','S03','S04','S05','S06'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
num = [];
diff=[];

%extract pupil size data
for kk = 2 : 2
    pupil_data_1 = [];
    pupil_data_2 = [];
    for i = 1 :1
        for j = 1 : 8
            main_folder = fullfile('F:\pupildata\Data_DI_wEYE\Data_DI_wEYE', subject{kk}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
    %        sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
            cd('MATs')
            cdt=dir([direction{j},'*_Dat_all.mat']);
            load(cdt.name)
            % we extract the time from trail_strat - 500 to trail_start +
            % 2000
            if j < 5 
                for kkk = 1 : 800
                    trial_st = tab(kkk,2) ;
                    aa=Dat_all(trial_st<=Dat_all(:,1) & Dat_all(:,1)<=trial_st+2100,:);
                    num = [num,size(aa,1)];
                    diff = [diff,max(aa(:,1))-min(aa(:,1))];
                    if size(aa,1) == 2101*samplingRateData/1000
                        [~,I]= sort(aa(:,1));
                        pupil_data_1 = [pupil_data_1;aa(I,4)'];
                    end                	
                end
            else
                for kkk = 1 : 800
                    trial_st = tab(kkk,2) ;
                    aa=Dat_all(trial_st<=Dat_all(:,1) & Dat_all(:,1)<=trial_st+2100,:);
%                     order = tab(i,2) < Dat_all(:,1) & tab(i,8) > Dat_all(:,1);
%                     trial= Dat_all(order,:);
                    
                    num = [num,size(aa,1)];
                    diff = [diff,max(aa(:,1))-min(aa(:,1))];
                    if size(aa,1) == 2101*samplingRateData/1000
                        [~,I]= sort(aa(:,1));
                        pupil_data_2 = [pupil_data_2;aa(I,4)'];
                    end                	
                end
            end
            sser1=sum(pupil_data_1 == 0 ,2);
            pupil_data_1=pupil_data_1(sser1==0,:);    
            sser2=sum(pupil_data_2 == 0 ,2);
            pupil_data_2=pupil_data_2(sser2==0,:);   
        end
    end
end

%%set baseline
baseline = [500,1000];
pupil_data_1no = preprocess_nor(baseline,pupil_data_1,samplingRateData);
pupil_data_2no = preprocess_nor(baseline,pupil_data_2,samplingRateData);

%% Visualize data
time = 0:2101*samplingRateData/1000-1;
figure(1)
plot(time,mean(pupil_data_1no,1),'r')
hold on
plot(time,mean(pupil_data_2no,1),'b')
xline(2600)
xline(3600)
title('Condition 1')
xlabel('time (ms)')  
ylabel('pupil size (% change from baseline')

figure(2)
plot(time,pupil_data_2no)
title('Condition 2')
xlabel('time (ms)')
ylabel('pupil size (% change from baseline')


%% define model
% Let's preallocate an empty model structure to represent a hypothetical
% task. We are going to use this model to generate artificial data for a
% single "subject".
taskmodel = pret_model();

% Let's assume that we have a sampling frequency of 1000 Hz and that the
% trials will be epoched from -500 to 3500 ms.
taskmodel.window = [0 2100.5];
taskmodel.samplerate = samplingRateData;
 
% Now let's suppose the trial sequence of this task is the following:
% 
%   A precue informs the observer that the trial has begun, the observer 
%   sees 2 stimuli 1000 and 1250 ms after the
%   precue, then a postcue 500 ms after the last stimulus instructs the
%   observer to respond. 
% 
% From this, we have 4 events; precue, stimulus 1, stimulus 2, postcue.
% Let's say the precue is time 0 ms, then these events occur at 0 ms, 1000
% ms, 1250 ms, and 1750 ms respectively.
taskmodel.eventtimes = [1300 1800];
taskmodel.eventlabels = { 'sti_on' 'sti_off'}; %optional

% We also have a response at the end of the trial. For the sake of 
% simplicity, let's say the observer always responds 1000 ms after the 
% postcue. Let's assume that there's a constant internal signal
% due to the cognitive workload associated with completing the task
% and making the decision. This constant internal signal would start at
% precue onset (0 ms) and last until time of response (2750 ms). In the
% context of this model, this would be a box regressor.
taskmodel.boxtimes = {[1300 2100]};
taskmodel.boxlabels = {'task'}; %optional

% Now we need to define the structure of the data using the model
% parameters. Let's say that the amplitude and latency of event-related
% pupil responses are expected to vary.
taskmodel.ampflag = true;
taskmodel.latflag = true;

% Let's also say we expect amplitude of the task-related pupil response and
% the tmax of the observer's pupil reponse to vary.
taskmodel.boxampflag = true;
taskmodel.tmaxflag = true;

% Now let's assume the baseline is perfectly stable (the y-intercept would
% always be 0) and that we don't expect any linear drift in pupil size
% during the trial (slope = 0).
taskmodel.yintflag = false;
taskmodel.slopeflag = false;

% Since the y-intercept and slope are always going to be 0, let's put that
% information into the structure.
taskmodel.yintval = 0;
taskmodel.slopeval = 0;

% Let's define boundaries for the parameters we will be fitting. Let's say
% event and box-related amplitudes will all be between 0 and 100 (percent 
% signal change from baseline), event latencies will be between -500 and 
% 500 ms, and tmax will be between 500 and 1500.
taskmodel.ampbounds = repmat([0;100],1,length(taskmodel.eventtimes));
taskmodel.latbounds = repmat([-500;500],1,length(taskmodel.eventtimes));
taskmodel.boxampbounds = [0;100];
taskmodel.tmaxbounds = [500;1500];

% So far, the model's specifications are common to all task conditions. Now
% let's define different models for each condition so that we can
% generate artificial data for each one. This is to simulate data that you
% might obtain in an experiment with three conditions.

%% organize data via pret_preprocess
% To fit pupil data you collected in an experiment, you would start here,
% with preprocessing. The data should be epoched by trial in a matrix of
% trials x time for each condition. pret_preprocess creates an sj
% (subject) structure with the data and metadata in a set format. If
% requested, it also performs baseline normalization and blink 
% interpolation on the trial data.

% The artificial data is already baseline normalized and has no blinks, so 
% we return the options structure and turn off those features
options = pret_preprocess();
options.normflag = false;
options.blinkflag = false;

% put the condition data matrices into a cell array
data = {pupil_data_1no pupil_data_2no};

% labels for each condition
condlabels = {'condition1' 'condition2'};

% other epoch info
samplerate = taskmodel.samplerate;
window = taskmodel.window;

sj = pret_preprocess(data,samplerate,window,condlabels,[],options);


%% create model for the data
% Pretending that we are naive to the model that we used to create our
% data, let's create a model to actually fit the data to.
model = pret_model();

% While the trial window of our task is from -500 to 3500 ms, here we are
% not interested in what's happening before 0. So
% let's set the model window to fit only to the region betweeen 0 and 3500
% ms (the cost function will only be evaluated along this interval).
model.window = [0 2100];

% We already know the sampling frequency.
model.samplerate = taskmodel.samplerate;

% We also know the event times of our task. Let's also say that we think 
% there will be a sustained internal signal from precue onset to response 
% time (0 to 2750 ms).
model.eventtimes = [1300 1800];
model.eventlabels = { 'sti_on' 'sti_off'}; %optional
model.boxtimes = {[1300 2100]};
model.boxlabels = {'task'}; %optional

% Let's say we want to fit a model with the following parameters: 
% event-related, amplitude, latency, task-related (box) amplitude, 
% and the tmax of the pupil response function. We turn the other parameters
% off.
model.yintflag = false;
model.slopeflag = false;

% Now let's define the bounds for the parameters we decided to fit. We do
% not have to give values for the y-intercept and slope because we are not
% fitting them.
model.ampbounds = repmat([0;100],1,length(model.eventtimes));
model.latbounds = repmat([-500;500],1,length(model.eventtimes));
model.boxampbounds = [0;100];
model.tmaxbounds = [500;1500];

% We need to fill in the values for the y-intercept and slope since we will
% not be fitting them as parameters.
model.yintval = 0;
model.slopeval = 0;

%% estimate model parameters via pret_estimate_sj
% Now let's perform the parameter estimation procedure on our subject data.
% The mean of each condition will be fit independently. For illustration, 
% let's run only 3 optimizations using one cpu worker (for more 
% information, see the help files of pret_estimate and pret_estimate_sj).
options = pret_estimate_sj(); 
options.pret_estimate.optimnum = 3;
% if you want to try fiting the parameters using single trials instead of the mean,
% use these lines (you'll want to turn off the optimization plots for this):
%   options.trialmode = 'single';
%   options.pret_estimate.pret_optim.optimplotflag = false;
wnum = 1;

sj = pret_estimate_sj(sj,model,wnum,options);

