function []=PF_Pilot_TF_Analysis(subj)
%% TF Analysis
% Input: Subject Number, generate plots (1/0)
fprintf('Starting Time Frequency Analysis Subj %i',subj)
singletrials=1; %(save single trials?)
artrej=0; % already remove artifact trials?

% Parameters (real analysis)
triggercodes={110;112;210;212}; % Warning Signals per condition (regular 800, catch 800, regular 850, catch 850)
timerange=[-300 1500];

paddingLength=500; % in ms
TF_wavFreqs=1:30; % Log Range
%TF_wavFreqs=2.^[0:1/6:5]; % Log Range

% Load Data
cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Preprocessed Data\PreprocessedEEG'
load(sprintf('EEG_PF_Pilot_Subj%i_pp.mat',subj))

if singletrials
        savefilename=sprintf('EEG_Pf_Pilot_Subj%i_TF_SingleTrials.mat',subj);
else
        savefilename=sprintf('EEG_SxA_Subj%i_TF_Results.mat',subj);
end


% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
srate=SDATA.info.sampling_rate;

% Segmentation
for c=1:size(triggercodes,1) % For conditions
    [segmentedData, isNotArtifact, timeVec]=segmentContEEGdata(triggercodes{c}, timerange+[-paddingLength paddingLength], data, triggers, artifacts, srate);

    % Check how many trials are artifact-free
    sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

    % Remove trials with artifacts
    if artrej==1
        segmentedData=segmentedData(:,:,isNotArtifact==1); 
    end

    % TF Analysis
    condResults_amp=[];
    condResults_phase=[];
    tic
    if ~singletrials
        parfor el=1:width(data)
            [wvlt_amp, wvlt_phase] = morletwave(TF_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate); % time points x frequnencies x trials
            inducedMat_a=squeeze(mean(wvlt_amp,3)); % average over trials
            condResults_amp(:,:,el)=inducedMat_a'; % reformat to: time points x frequencies x electrodes
            inducedMat_p=squeeze(mean(wvlt_phase,3)); % average over trials
            condResults_phase(:,:,el)=inducedMat_p'; % reformat to: time points x frequencies x electrodes
        end
    else
        parfor el=1:width(data)
            [wvlt_amp, wvlt_phase] = morletwave(TF_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate); % time points x frequnencies x trials
            inducedMat_a=wvlt_amp; % do not average over trials
            condResults_amp(:,:,:,el)=permute(inducedMat_a,[2,1,3]); % reformat to: time points x frequencies x trials x electrodes
            inducedMat_p=wvlt_phase; % do not average over trials
            condResults_phase(:,:,:,el)=permute(inducedMat_p,[2,1,3]); % reformat to: time points x frequencies x trials x electrodes
        end
    end
    toc

    % Remove Padding
    if ~singletrials
        condResults_amp=condResults_amp(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:);
        condResults_phase=condResults_phase(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:);
    else
        condResults_amp=condResults_amp(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:,:);
        condResults_phase=condResults_phase(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:,:);
    end

    timeVec=timeVec(timeVec>=timerange(1) & timeVec<=timerange(2));
    TF_Results_amp{c}=condResults_amp; % Time Points, Frequencies, Trials, Electrodes
    TF_Results_phase{c}=condResults_phase; % Time Points, Frequencies, Trials, Electrodes
    TF_timeVecTotal{c}=timeVec;
    TF_NotArtifact{c}=isNotArtifact;
    clear condResults segmentedData timeVec
end
cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\TimeFrequency'
if ~singletrials
    save(savefilename, 'TF_Results_amp','TF_Results_phase', 'TF_timeVecTotal', 'TF_wavFreqs','TF_NotArtifact');
else
    TF_Results_Trial_amp=TF_Results_amp;
    TF_Results_Trial_phase=TF_Results_phase;
    TF_trial_timeVec=TF_timeVecTotal;
    save(savefilename, 'TF_Results_Trial_amp','TF_Results_Trial_phase','TF_trial_timeVec','TF_NotArtifact','-v7.3');
end
fprintf('TF Results Saved - Subj %i',subj)
end