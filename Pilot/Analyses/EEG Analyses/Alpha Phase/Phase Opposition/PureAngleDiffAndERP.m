%% Plot Alpha ERP Time Course and Angle Difference
clear
clc

% What to plot?
subj=[1:9];
alpha_range=[8 12]; % upper and lower bound of frequency range
gl_plots=1; % also plot average across all participants in subj?
% Params
elec=[25:30 62:64]; % occipital electrodes

% Segmentation Parameters
trigcodes_target={110;210}; % trigger codes for target trials (at regular time points) - currently cue-codes for each condition
trigcodes_catch={112;212}; % trigger codes for catch trials only (cue-codes) for each condition
timerange=[-200 1100]; % 0 is cue
artrej=1; % reject artifact trials?

% Preallocate
GL_cond1_target=0;
GL_cond2_target=0;
GL_cond1_catch=0;
GL_cond2_catch=0;

% Load phase data and calculate angle distance
for s=1:length(subj)
    cd  'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Alpha Phase'
    for t=1:2 % for the two target conditions (target and catch)
        nexttile;
        if t==1
            % Load Target Trials
            loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Target.mat',subj(s));
            load(loadfilename1, "Alpha_SingleTrials","ITPCalpha_timevec")
        elseif t==2
            % Load Catch Trials
            loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Catch.mat',subj(s));
            load(loadfilename1, "Alpha_SingleTrials","ITPCalpha_timevec")
        end
        % Prep (separate by condition, select electrodes and average)
        cond1_angles=Alpha_SingleTrials{1,1};
        cond2_angles=Alpha_SingleTrials{1,2};
        cond1_angles=squeeze(circ_mean(cond1_angles(:,elec,:),[],2)); % average across electrodes
        cond2_angles=squeeze(circ_mean(cond2_angles(:,elec,:),[],2));
        cond1_mean=circ_mean(cond1_angles,[],2); % average across trials
        cond2_mean=circ_mean(cond2_angles,[],2);
        anglediff(s,t,:)=circ_dist(cond1_mean,cond2_mean);
    end
    % Get ERP Data
    loadfilename=sprintf('EEG_PF_Pilot_Subj%i_pp.mat',subj(s));
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Preprocessed Data\PreprocessedEEG'
    load(loadfilename)
    all_data=SDATA.data;
    artifacts=SDATA.metadata.artifacts;
    triggers=SDATA.events.triggerChannel;
    srate=SDATA.info.sampling_rate;
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

for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(1:size(anglediff,3),squeeze(anglediff(s,1,:)))
    hold on
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_target(s,1,:)))
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_target(s,2,:)))
    nexttile;
    title("Catch Trials")
    plot(1:size(anglediff,3),squeeze(anglediff(s,2,:)))
    hold on
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_catch(s,1,:)))
    plot(1:size(Alpha_ERP_catch,3),squeeze(Alpha_ERP_catch(s,2,:)))
end
%% Plot
for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(1:size(anglediff,3),squeeze(anglediff(s,1,:)))
    hold on
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_target(s,1,:)))
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_target(s,2,:)))
    xline(0,"WS")
    xline(800,"Target")
    nexttile;
    title("Catch Trials")
    plot(1:size(anglediff,3),squeeze(anglediff(s,2,:)))
    hold on
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_catch(s,1,:)))
    plot(1:size(Alpha_ERP_catch,3),squeeze(Alpha_ERP_catch(s,2,:)))
    xline(0,"WS")
    xline(800,"Target")
end
xline(0)
for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(1:size(anglediff,3),squeeze(anglediff(s,1,:)))
    hold on
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_target(s,1,:)))
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_target(s,2,:)))
    xline(0)
    xline(800)
    nexttile;
    title("Catch Trials")
    plot(1:size(anglediff,3),squeeze(anglediff(s,2,:)))
    hold on
    plot(1:size(Alpha_ERP_target,3),squeeze(Alpha_ERP_catch(s,1,:)))
    plot(1:size(Alpha_ERP_catch,3),squeeze(Alpha_ERP_catch(s,2,:)))
    xline(0)
    xline(800)
end
for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(timeVec,squeeze(anglediff(s,1,:)))
    hold on
    plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)))
    plot(timeVec,squeeze(Alpha_ERP_target(s,2,:)))
    xline(0)
    xline(800)
    nexttile;
    title("Catch Trials")
    plot(timeVec,squeeze(anglediff(s,2,:)))
    hold on
    plot(timeVec,squeeze(Alpha_ERP_catch(s,1,:)))
    plot(timeVec,squeeze(Alpha_ERP_catch(s,2,:)))
    xline(0)
    xline(800)
end
for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(timeVec,squeeze(anglediff(s,1,:)),'LineWidth',2)
    hold on
    plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
    plot(timeVec,squeeze(Alpha_ERP_target(s,2,:)),'LineWidth',2)
    xline(0)
    xline(800)
    nexttile;
    title("Catch Trials")
    plot(timeVec,squeeze(anglediff(s,2,:)),'LineWidth',2)
    hold on
    plot(timeVec,squeeze(Alpha_ERP_catch(s,1,:)),'LineWidth',2)
    plot(timeVec,squeeze(Alpha_ERP_catch(s,2,:)),'LineWidth',2)
    xline(0)
    xline(800)
end
title("Target Trials")
plot(timeVec,squeeze(anglediff(s,1,:)),'LineWidth',2)
hold on
plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
legend(["Angle Difference","ERP 800","ERP 850"])
legend({"Angle Difference","ERP 800","ERP 850"})
legend("Angle Difference","ERP 800","ERP 850")
nexttile;
title("Target Trials")
plot(timeVec,squeeze(anglediff(s,1,:)),'LineWidth',2)
hold on
plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
plot(timeVec,squeeze(Alpha_ERP_target(s,2,:)),'LineWidth',2)
xline(0)
xline(800)
legend("Angle Difference","ERP 800","ERP 850")
plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
nexttile;
title("Target Trials")
plot(timeVec,squeeze(anglediff(s,1,:)),'LineWidth',2)
hold on
plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
plot(timeVec,squeeze(Alpha_ERP_target(s,2,:)),'LineWidth',2)
xline(0)
xline(800)
legend("Angle Difference","ERP 800","ERP 850")
for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(timeVec,squeeze(anglediff(s,1,:)),'LineWidth',2)
    hold on
    plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
    plot(timeVec,squeeze(Alpha_ERP_target(s,2,:)),'LineWidth',2)
    xline(0)
    xline(800)
    legend("Angle Difference","ERP 800","ERP 850")
    nexttile;
    title("Catch Trials")
    plot(timeVec,squeeze(anglediff(s,2,:)),'LineWidth',2)
    hold on
    plot(timeVec,squeeze(Alpha_ERP_catch(s,1,:)),'LineWidth',2)
    plot(timeVec,squeeze(Alpha_ERP_catch(s,2,:)),'LineWidth',2)
    xline(0)
    xline(800)
    legend("Angle Difference","ERP 800","ERP 850")
end

for s=1:length(subj)
    figure; t=tiledlayout('flow');
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    nexttile;
    title("Target Trials")
    plot(timeVec,squeeze(anglediff(s,1,:)),'LineWidth',2)
    hold on
    plot(timeVec,squeeze(Alpha_ERP_target(s,1,:)),'LineWidth',2)
    plot(timeVec,squeeze(Alpha_ERP_target(s,2,:)),'LineWidth',2)
    xline(0)
    xline(800)
    xline(850)
    legend("Angle Difference","ERP 800","ERP 850")
    nexttile;
    title("Catch Trials")
    plot(timeVec,squeeze(anglediff(s,2,:)),'LineWidth',2)
    hold on
    plot(timeVec,squeeze(Alpha_ERP_catch(s,1,:)),'LineWidth',2)
    plot(timeVec,squeeze(Alpha_ERP_catch(s,2,:)),'LineWidth',2)
    xline(0)
    xline(800)
    legend("Angle Difference","ERP 800","ERP 850")
end