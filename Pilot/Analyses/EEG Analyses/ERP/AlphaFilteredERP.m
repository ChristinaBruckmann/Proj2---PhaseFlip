%% Alpha Band Passed ERP
clear
clc

subj=[9];
alpha_range=[8 12]; % upper and lower bound of frequency range
elec=[25:30 62:64]; % occipital electrodes
noseref=0; % rereference all the data to the nose?

% Segmentation Parameters
trigcodes_target={110;210}; % trigger codes for target trials (at regular time points) - currently cue-codes for each condition
trigcodes_catch={112;212}; % trigger codes for catch trials only (cue-codes) for each condition
timerange=[-200 1100]; % 0 is cue
artrej=1; % reject artifact trials?

for s=1:length(subj)
    % Load Data
    loadfilename=sprintf('EEG_PF_Pilot_Subj%i_pp.mat',subj(s));
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Preprocessed Data\PreprocessedEEG'
    load(loadfilename)

    all_data=SDATA.data;
    artifacts=SDATA.metadata.artifacts;
    triggers=SDATA.events.triggerChannel;
    srate=SDATA.info.sampling_rate;

    % Re-reference to the nose?
    if noseref

    end

    % Filter
    [alpha_whole_ts] = bandPassFilter(alpha_range(1),alpha_range(2),all_data,srate); % filter around alpha range

    % Segment
    for tc=1:length(trigcodes_target)
        [alpha_st_target, isNotArtifact, timeVec]=segmentContEEGdata(trigcodes_target{tc} , timerange,...
            alpha_whole_ts , triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        if artrej==1
            Alpha_SingleTrials_target=alpha_st_target(:,:,logical(isNotArtifact));
            Alpha_ERP_target(s,tc,:,:)=mean(Alpha_SingleTrials_target,3);
        end
    end

    for tc=1:length(trigcodes_catch)
        [alpha_st_catch, isNotArtifact, timeVec]=segmentContEEGdata(trigcodes_catch{tc} , timerange,...
            alpha_whole_ts , triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        if artrej==1
            Alpha_SingleTrials_catch=alpha_st_catch(:,:,logical(isNotArtifact));
            Alpha_ERP_catch(s,tc,:,:)=mean(Alpha_SingleTrials_catch,3);
        end
    end
end
% Select and average across electrodes
Alpha_ERP_target=mean(Alpha_ERP_target(:,:,:,elec),4);
Alpha_ERP_catch=mean(Alpha_ERP_catch(:,:,:,elec),4);

% Average across subjects
Alpha_ERP_target=squeeze(mean(Alpha_ERP_target,1));
Alpha_ERP_catch=squeeze(mean(Alpha_ERP_catch,1));

% Aggregate for plotting
alpha_data(1,:,:)=Alpha_ERP_target;
alpha_data(2,:,:)=Alpha_ERP_catch;

% Plot
colourvec=[0.4940 0.1840 0.5560 1;
    0.4660 0.6740 0.1880 1];


maintitle=sprintf("Subject %i",subj);
figure;tiledlayout('flow');title('Valid Target')
for  t=1:2 % two target types (target and catch)
    nexttile;
    for c=1:2 % two timing conditions (800 and 850)
        plot(timeVec,squeeze(alpha_data(t,c,:)),"LineWidth",4,"Color",colourvec(c,:))
        hold on;
        if t==1
            title('Valid Target')
        else
            title('Catch Trials')
        end
    end
    xline(0)
    xline(800)
    xline(850)
    legend("800","850")
    box('off')
end

