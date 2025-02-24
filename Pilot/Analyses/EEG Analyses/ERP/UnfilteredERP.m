%% Unfiltered ERP
clear
clc

subj=[7];
elec=[25:30 62:64]; % occipital electrodes

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

    % Segment
    for tc=1:length(trigcodes_target)
        [alpha_st_target, isNotArtifact, timeVec]=segmentContEEGdata(trigcodes_target{tc} , timerange,...
            all_data , triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        if artrej==1
            SingleTrials_target=alpha_st_target(:,:,logical(isNotArtifact));
            ERP_target(s,tc,:,:,:)=mean(SingleTrials_target,3);
        end
    end

    for tc=1:length(trigcodes_catch)
        [alpha_st_catch, isNotArtifact, timeVec]=segmentContEEGdata(trigcodes_catch{tc} , timerange,...
            all_data, triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        if artrej==1
            SingleTrials_catch=alpha_st_catch(:,:,logical(isNotArtifact));
            ERP_catch(s,tc,:,:,:)=mean(SingleTrials_catch,3);
        end
    end
end

% Select and average across electrodes
ERP_target=mean(ERP_target(:,:,:,elec),4);
ERP_catch=mean(ERP_catch(:,:,:,elec),4);

% Average across subjects
ERP_target=squeeze(mean(ERP_target,1));
ERP_catch=squeeze(mean(ERP_catch,1));

% Lowpassfilter
lpf_cutoff=20;
for tc=1:2
ERP_target(tc,:)=LPF(ERP_target(tc,:), srate, lpf_cutoff);
ERP_catch(tc,:)=LPF(ERP_catch(tc,:), srate, lpf_cutoff);
end

% Aggregate for plotting
final_data(1,:,:)=ERP_target;
final_data(2,:,:)=ERP_catch;

% Plot
colourvec=[0.4940 0.1840 0.5560 1;
    0.4660 0.6740 0.1880 1];

figure;tiledlayout('flow');
for  t=1:2 % two target types (target and catch)
    nexttile;
    for c=1:2 % two timing conditions (800 and 850)
        plot(timeVec,squeeze(final_data(t,c,:)),"LineWidth",4,"Color",colourvec(c,:))
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
