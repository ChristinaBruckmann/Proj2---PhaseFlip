%% Calculates the Bifurcation Index for each frequency and time point between the 800 and 850 conditions.
% According to van Rullen 2016 and Busch 2009
% PBI = (ITCa − ITC_all)(ITCb − ITC_all)
clear
clc

% Settings
subj=[1:9];
freq_range=1:30;
elec=[25:30 62:64]; % occipital electrodes
trialtype=[1 2]; % 1 - target trials, 2 - catch trials, 1 2 plot both separately

% Optional Cleaning
after_reminder=0; % only trials immediately after reminder interval
threshold_only=1; % only trials with a contrast around the participants threshold (irrelevant for catch)
after_learning=1; % only trials after the initial 10 trials of main experiment (first 20 practice trials are excluded anyway)
correct_only=1; % only trials with correct responses (irrelevant for catch)

for s=1:length(subj)
    % Load Data
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\TimeFrequency'
    load(sprintf('EEG_Pf_Pilot_Subj%i_TF_SingleTrials.mat',subj(s)),'TF_Results_Trial_phase', 'TF_trial_timeVec', 'TF_NotArtifact');
     
    % Get time vector (same for all conditions)
    timeVec=TF_trial_timeVec{1,1};
    clear TF_trial_timeVec

    % Load Behavioural Data
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Raw Data\Raw Behaviour'
    load(sprintf('Pilot_PhaseFlip_Subj%i.mat',s))

    % Determine trial type and select data
    for c=1:2 % conditions (1- 800 2-850)
        for t=1:length(trialtype)
            % Get relevant behavioural data for trial rejection
            data=subresults.data(subresults.data{:,"Condition"}==c,:); % choose correct condition

            % Optional Cleaning:
            main_include=ones(height(data),1);

            % Create index to select only trials after reminders
            if after_reminder
                after_reminder_idx=zeros(height(data),1); % preallocate logical index
                after_reminder_idx(1:4:height(data))=1; % reminder every 4 trials, mark as logical 1
                main_include=main_include & after_reminder_idx; % add to main logical index
            end

            % Create index to select only trials where the contrast was around the threshold
            if threshold_only && trialtype(t)==1 % threshold does not apply to catch trials
                threshold_std=std(data{10:end,"Gabor Strength"}); % exclude first 10 because of rough calibration in the beginning
                threshold_mean=mean(data{10:end,"Gabor Strength"}); % exclude first 10 because of rough calibration in the beginning
                threshold_range=[(threshold_mean-2*threshold_std) (threshold_mean+2*threshold_std)]; % get range +- 2 stds from mean
                within_thresh_bounds_idx=data{:,"Gabor Strength"}>min(threshold_range) & data{:,"Gabor Strength"}<max(threshold_range);
                main_include=main_include & within_thresh_bounds_idx;
            end

            % Exclude the first trials for learning
            if after_learning
                after_learning_idx=ones(height(data),1); % preallocate logical index
                after_learning_idx(1:10)=0; % preallocate logical index
                main_include=main_include & after_learning_idx;
            end

            % Create index to only include correct trials
            if correct_only && trialtype(t)==1 % correct does not apply to catch trials
                correct_response_idx=zeros(height(data),1); % preallocate logical index
                correct_response_idx(data{:,"Correct/Incorrect"}==1,:)=1;
                main_include=main_include & correct_response_idx;
            end

            % Select only the actual trial type (condition and catch/target) from the selection vector
            trialtype_idx=zeros(height(data),1); % preallocate logical index
            if trialtype(t)==1
                trialtype_idx(data{:,"Trial Type"}==1,:)=1; % only select target trials
            elseif trialtype(t)==2
                trialtype_idx(data{:,"Trial Type"}==5,:)=1; % only select target trials
            else
                error('unknown trial type')
            end
            main_include=main_include(logical(trialtype_idx));

            % Trials without artifacts index
            if c==1
                artifact_free_idx=TF_NotArtifact{1,trialtype(t)};
            else
                artifact_free_idx=TF_NotArtifact{1,trialtype(t)+2};
            end

            if trialtype(t)==1
                artifact_free_idx=artifact_free_idx(21:end); %target trials have 20 EEG trials that are practice that are not recorded in the behavioural data, so these are being cut here
            end

            main_include=main_include & artifact_free_idx;

            % Select EEG Data
            if c==1
                currEEG_data=TF_Results_Trial_phase{1, trialtype(t)};
            elseif c==2
                currEEG_data=TF_Results_Trial_phase{1, trialtype(t)+c};
            end

            % Filter EEG Data
            if trialtype(t)==1 %target trials have 20 EEG trials that are practice that are not recorded in the behavioural data, so these are being cut here
                currEEG_data=currEEG_data(:,:,21:end,:);
            end

            clean_data=currEEG_data(:,freq_range,:,elec); % select occipital electrodes and alpha freqs only
            clean_data=clean_data(:,:,main_include,:); % choose selected trials only (based on selection criteria above)

            % Average Data
            mean_data=squeeze(mean(clean_data,4)); % average across electrodes

            % Save Data
            if c==1
                GL_TF_Phase{s,trialtype(t)}=mean_data;
            elseif c==2
                GL_TF_Phase{s, trialtype(t)+c}=mean_data;
            end
        end
    end
end

% Save GL File
cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Bifurcation'
save('GL_TF_Phase','GL_TF_Phase','timeVec','after_reminder','threshold_only','correct_only','after_reminder')


%% Calculate overall ITC (sum across trials: (angle / amplitude) then divide by number of trials)

for s=1:size(GL_TF_Phase,1) % For each subject

    f=figure; tiledlayout('flow');
    %title(f,sprintf('Subject %i',s))

    for t=1:(size(GL_TF_Phase,2)/2) % for trialtype (target and catch)

        % Select Data for subj x condition x trialtype
        data800=GL_TF_Phase{s,t};
        data850=GL_TF_Phase{s, t+2};

        % Calculate Phase Opposition Values
        [p_circWW, p_POS, p_zPOS]=PhaseOpposition(data800, data850);

        % Plot
         nexttile;
         imagesc(p_circWW');
         colorbar;

         nexttile;
         imagesc(p_POS');
         colorbar;

         nexttile;
         imagesc(p_zPOS');
         colorbar;

        % Save
        subj_res(1,:,:)=p_circWW;
        subj_res(2,:,:)=p_POS;
        subj_res(3,:,:)=p_zPOS;
        GL_PhaseOppositionResults(s,t,:,:,:)=subj_res; % subject x trial type (target or catch) x type of analysis (circWW,POS,zPOS) x time points x frequencies
    end
end

% Save
cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Bifurcation'
save('GL_PhaseOppositionIdx','GL_PhaseOppositionResults','timeVec','after_reminder','threshold_only','correct_only','after_reminder')

%% Group Level Plot
% Average across subjects
GL_POR_mean=squeeze(mean(GL_PhaseOppositionResults,1));  % trial type x analysis x time points x frequencies
f=figure; tiledlayout('flow');

for t=1:size(GL_POR_mean,1) % for trialtype (target and catch)
    % Plot
    nexttile;
    imagesc(timeVec, 1:30, squeeze(GL_POR_mean(t,1,:,:))'); axis xy; 
    colorbar;
    xline(0)
    xline(800)
    xline(850)
    if t==1
        title("Target Trials - circ WW")
    else
        title("Catch Trials - circ WW")
    end

    nexttile;
    imagesc(timeVec, 1:30, squeeze(GL_POR_mean(t,2,:,:))'); axis xy; 
    colorbar;
    xline(0)
    xline(800)
    xline(850)
     if t==1
        title("Target Trials - Phase Opposition Sum")
    else
        title("Catch Trials - Phase Opposition Sum")
    end

    nexttile;
    imagesc(timeVec, 1:30, squeeze(GL_POR_mean(t,3,:,:))'); axis xy; 
    colorbar;
    xline(0)
    xline(800)
    xline(850)
     if t==1
        title("Target Trials - z-scored POS")
    else
        title("Catch Trials - z-scored POS")
    end
end
