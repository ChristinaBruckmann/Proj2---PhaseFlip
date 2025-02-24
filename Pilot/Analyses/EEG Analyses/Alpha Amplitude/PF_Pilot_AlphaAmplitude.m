%% Plot Alpha Amplitude
% Requires TF results (single trials, no artifact rejected) for every participant

subj=[1:9];
alpha_range=8:12; % upper and lower bound of frequency range
elec=[25:30 62:64]; % occipital electrodes
trialtype=[1 2]; % 1 - target trials, 2 - catch trials, 1 2 plot both separately

% Optional Cleaning
after_reminder=0; % only trials immediately after reminder interval
threshold_only=0; % only trials with a contrast around the participants threshold (irrelevant for catch)
after_learning=0; % only trials after the initial 10 trials of main experiment (first 20 practice trials are excluded anyway)
correct_only=0; % only trials with correct responses (irrelevant for catch)

for s=1:length(subj)
    % Load Data
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\TimeFrequency'
    load(sprintf('EEG_Pf_Pilot_Subj%i_TF_SingleTrials.mat',subj(s)))

    timeVec=TF_trial_timeVec{1,1}; % simplify time vec, same for all trials anyway
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
                currEEG_data=TF_Results_Trial_amp{1, trialtype(t)};
            elseif c==2
                currEEG_data=TF_Results_Trial_amp{1, trialtype(t)+c};
            end

            % Filter EEG Data
            if trialtype(t)==1 %target trials have 20 EEG trials that are practice that are not recorded in the behavioural data, so these are being cut here
                currEEG_data=currEEG_data(:,:,21:end,:);
            end

            clean_data=currEEG_data(:,alpha_range,:,elec); % select occipital electrodes and alpha freqs only
            clean_data=clean_data(:,:,main_include,:); % choose selected trials only (based on selection criteria above)

            % Average Data
            mean_data=squeeze(mean(clean_data,2)); % average across frequencies
            mean_data=squeeze(mean(mean_data,2)); % average across trials
            mean_data=squeeze(mean(mean_data,2)); % average across electrodes

            % Save Data
            GL_AlphaAmp(subj(s),t,c,:)=mean_data;
        end
    end
end

%% Plot Group Level Results
bc=1; % baseline correct?
bc_tw=[0 100]; % time window for baseline correction (in relation to WS)

% Average across participants
GL_AlphaAmp_mean=squeeze(mean(GL_AlphaAmp,1));

% Baseline Correct
if bc
    bc_idx=timeVec>=bc_tw(1)&timeVec<=bc_tw(2); % index the baseline time window
    baseline=mean(GL_AlphaAmp_mean(:,:,bc_idx),3); % calculate mean across that time window
    GL_AlphaAmp_mean=GL_AlphaAmp_mean-baseline; % subtract mean from data
end

% Plot
figure;
for c=1:2
    for t=1:length(trialtype)

        % Assign line colour (800 blue, 850 red, target dark, catch light)
        if c==1
            if trialtype(t)==1
                colorvec=[0 0.4470 0.7410]; % dark blue
            else
                colorvec=[0.3010 0.7450 0.9330]; % light blue
            end
        else
            if trialtype(t)==2
                colorvec=[0.6350 0.0780 0.1840]; % red
            else
                colorvec= [0.8500 0.3250 0.0980]; % orange
            end
        end

        plot(timeVec,squeeze(GL_AlphaAmp_mean(t,c,:)),'Color',colorvec,'LineWidth',3)
        hold on
    end
end

% Make plot nice
xline(800,'Color',[0 0.4470 0.7410],'LineWidth',2) % 800ms target
xline(850,'Color',[0.6350 0.0780 0.1840],'LineWidth',2) % 850 ms target
xline(0,'Label','Warning Signal')
legend("800 Target", "800 Catch","850 Target","850 Catch")
title("Pre-Target Suppression of Alpha Amplitude")

%% Jackknift Plots
figure;tiledlayout('flow')

bc=1; % baseline correct?
bc_tw=[0 100]; % time window for baseline correction (in relation to WS)

for s=1:size(GL_AlphaAmp,1)
    % Remove current subject
    jackknife_data=GL_AlphaAmp;
    jackknife_data(s,:,:,:)=[];
    % Calculate mean over all other subjects
    jackknife_mean =squeeze(mean( jackknife_data,1));

    % Baseline Correct
    if bc
        bc_idx=timeVec>=bc_tw(1)&timeVec<=bc_tw(2); % index the baseline time window
        baseline=mean(jackknife_mean(:,:,bc_idx),3); % calculate mean across that time window
        jackknife_mean=jackknife_mean-baseline; % subtract mean from data
    end

    % Plot
    nexttile;
    for c=1:2
        for t=1:length(trialtype)

            % Assign line colour (800 blue, 850 red, target dark, catch light)
            if c==1
                if trialtype(t)==1
                    colorvec=[0 0.4470 0.7410]; % dark blue
                else
                    colorvec=[0.3010 0.7450 0.9330]; % light blue
                end
            else
                if trialtype(t)==2
                    colorvec=[0.6350 0.0780 0.1840]; % red
                else
                    colorvec= [0.8500 0.3250 0.0980]; % orange
                end
            end

            plot(timeVec,squeeze(jackknife_mean(t,c,:)),'Color',colorvec,'LineWidth',3)
            hold on
        end
    end

    % Make plot nice
    xline(800,'Color',[0 0.4470 0.7410],'LineWidth',2) % 800ms target
    xline(850,'Color',[0.6350 0.0780 0.1840],'LineWidth',2) % 850 ms target
    xline(0,'Label','Warning Signal')
    legend("800 Target", "800 Catch","850 Target","850 Catch")
    title(sprintf("Without Subj %i",s))
end