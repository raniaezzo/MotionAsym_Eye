




function [temporal_idx1, latency, RT_output, accuracy_output, trialmatOutput, dataCheck, supTimes, postOffsetTimes, selectedTrials] = ...
    compute_latency(data, fieldname, suppressionTimes, postsupPeakTimes, rtData, trialData, stimonset)
    
    if ~isstruct(data)
        error('data must be a struct.')
    else
        if isempty(fieldname)
            names = fieldnames(data);
        else
            names = {fieldname};
        end
    end

    %[rows, cols] = size(matrix);
    %filtered_matrix = zeros(rows, cols);

    trialmatOutput = [];

    for fi=1:length(names)

        temporal_idx1 = []; % I moved these into the loop
        temporal_idx2 = [];
        latency = [];
    
        RT_output = [];
        accuracy_output = []; %%
        supTimes = [];
        postOffsetTimes = [];

        matrix = data.(names{fi});
        [rws, ~] = size(matrix);
        filtered_matrix = zeros(size(matrix));

        rtMat = rtData.(names{fi});
        trialMat = trialData.(names{fi});

        selectedTrials = []; % which rows "qualify"

        % find the post stim inflection point (if not already given)
        window_size = 250; % Window size for the moving average
        smoothed_data = smoothdata(nanmean(matrix,1), 'movmean', window_size);
        
        decreasing_index = find(diff(smoothed_data(1800:end)) < 0, 1);
        postOffsetDip = decreasing_index+1800;

        if (length(suppressionTimes) == 1)
            suppressionTime = ones(rws,1)*suppressionTimes;
            postsupPeakTime = ones(rws,1)*postOffsetDip; % this is the mean peak for ALL data
        else
            temp = mean(suppressionTimes);
            %suppressionTime = ones(rws,1)*temp; % just take the mean?
            suppressionTime=suppressionTimes;
            postsupPeakTime = postsupPeakTimes;
        end

        

%         figure
%         plot(smoothed_data)
%         hold on
%         xline(1300)
%         hold on
%         xline(1800)
%         hold on
%         xline(decreasing_index+1800, 'b--')

        for i = 1:rws

%             if ~ismember(trialMat(i,10),[3,7])
%                 continue
%             end

            leniency = suppressionTime(i); % allow pre stim onset to leak in

            row = matrix(i,:);
   
            idx1 = []; idx2 = [];
            try
                % first MS after STIM ONSET
                %idx1 = find(row(stimonset+leniency:end) == 1, 1);
                idx1 = find(row(stimonset+leniency:floor(rtMat(i)+1300)) == 1, 1); %postsupPeakTime(i))==1,1, 'first'); % 
                %%stimonset+leniency+1000)== 1, 1); %floor(rtMat(i))) == 1, 1); % +500) == 1, 1); %postOffsetDip)==1,1); %
                idx2 = find(row(stimonset-1000:stimonset+leniency) == 1, 1, 'last'); %'last');
            catch
                idx1 = find(row(stimonset+leniency:end) == 1, 1);
                %idx2 = find(row(1:stimonset+leniency) == 1, 1);
                idx2 = find(row(1:stimonset+leniency) == 1, 1, 'last'); %'last');
            end
            
    
            if ~isempty(idx1) && ~isnan(rtMat(i)) %&& ~isempty(idx2)
%                 %filtered_matrix(i, stimonset+leniency + idx1) = 2; % post stim onset
%                 %filtered_matrix(i, idx2) = 1; % pre stim onset
                temporal_idx1 = [temporal_idx1, stimonset+leniency + idx1];
                temporal_idx2 = [temporal_idx2, idx2];
                latency = [latency, (stimonset+leniency + idx1) - idx2];
                selectedTrials = [selectedTrials i];
                supTimes = [supTimes; stimonset+leniency];
                postOffsetTimes = [postOffsetTimes; postsupPeakTime(i)];
            else
                % do nothing?
            end
        end

        RT_output = [RT_output; rtMat(selectedTrials)];
        accuracy_output = [accuracy_output; trialMat(selectedTrials, 14)];
        trialmatOutput{fi} = trialMat(selectedTrials,:);
        dataCheck = data.allMSData(selectedTrials,:);
%    
%         easyTrials = ismember(trialmatOutput{fi}(:,10), 5:8);
%         hardTrials = ismember(trialmatOutput{fi}(:,10), 1:4);
%         rate_easy = nansum(dataCheck(easyTrials,:),1);
%         rate_hard = nansum(dataCheck(hardTrials,:),1);
%         figure
%         plot(rate_easy,'g')
%         hold on
%         plot(rate_hard, 'r')
%         hold on
%         xline(leniency+stimonset, 'k')
%         hold on
%         xline(postOffsetDip, 'k')
%         hold on
%         xline(1300, 'b')
%         hold on
%         xline(1800, 'b')

    end
end
