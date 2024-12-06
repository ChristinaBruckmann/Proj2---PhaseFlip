%% Single Subject EEG Alpha Phase Analysis (Phase Flip)
% Basic Alpha Analysis, provides the base for all other alpha codes of this project. Needs to be run for each participant before other analyses can be done.
clear
clc

subj=3; % vector with multiple subjects also works
noisy_nose=[]; % subjects from subj with a noisy nose (which should be referenced to mastoids instead)

% Segmentation Parameters
triggercodes_target={110;210}; % trigger codes for target trials (at regular time points) - currently cue-codes for each condition
triggercodes_catch={112;212}; % trigger codes for catch trials only (cue-codes) for each condition
timerange=[-200 1100]; % 0 is cue

% % Get Alpha Phase for the entire time series
savefilename1='EEG_PF_Pilot_Subj%i_AlphaPhaseWholeTS.mat';
PF_AlphaPhaseExtraction(subj,noisy_nose,savefilename1)

% Cut alpha phase around triggers
savefilename2='EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Target.mat';
savefilename3='EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Catch.mat';
PF_AlphaPhaseSegmentation(subj,triggercodes_target,timerange,savefilename2)
PF_AlphaPhaseSegmentation(subj,triggercodes_catch,timerange,savefilename3)

%% Extract Alpha Phase from the entire time series
% Individual Trial data gets saved to that they can all be plotted and compared.

function PF_AlphaPhaseExtraction(subj, noisy_ref,savefn)

for s=1:length(subj)
    disp('Starting Alpha Phase Analysis')

    % Load Data
    loadfilename=sprintf('EEG_PF_Pilot_Subj%i_pp.mat',subj(s));
    savefilename=sprintf(savefn,subj(s));
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Preprocessed Data\PreprocessedEEG'
    load(loadfilename)

    % Parameters
    ITPCalpha_wavFreqs=[10*2^-0.3 10*2^0.3]; % confirm that this makes sense
    srate=SDATA.info.sampling_rate;
    ITPCalpha_elec= 1:64;

    % Re-Reference for different clusters
    if ismember(subj(s),noisy_ref) % has a noisy nose, always refernce to mastoids
        data_occ=referenceContEEGdata(SDATA.data,[69 71]);
        fprintf("Subject %i referenced to mastoids.",subj(s))
    else
        data_occ=referenceContEEGdata(SDATA.data,71); %nose
    end

    % band-pass filter 0.5-2 Hz butterworth, 24dB/octave
    bpFilteredData_occ = bandPassFilter(min(ITPCalpha_wavFreqs),max(ITPCalpha_wavFreqs),data_occ(:,ITPCalpha_elec),srate);

    % Hilbert Transform
    hil_occ = hilbert(bpFilteredData_occ);

    % Extract Phase
    alpha_phase_whole_ts = angle(hil_occ); % phase

    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Alpha Phase'
    save(savefilename, 'alpha_phase_whole_ts','-v7.3')
    disp('Saved Alpha Phase Analysis')
end
end

%% Alpha Segmentation into Trials
% Segmentation of the whole time-series alpha values into trials (or around selected triggers)
function PF_AlphaPhaseSegmentation(subj,trigcodes,time,savefn)


for s=1:length(subj)

    loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseWholeTS.mat',subj(s));
    loadfilename2=sprintf('EEG_PF_Pilot_Subj%i_pp.mat',subj(s));
    savefilename=sprintf(savefn,subj(s));

    % Load
    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Alpha Phase'
    load(loadfilename1,'alpha_phase_whole_ts') % Load alpha phase

    cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Preprocessed Data\PreprocessedEEG'
    load(loadfilename2,'SDATA')

    artifacts=SDATA.metadata.artifacts;
    triggers=SDATA.events.triggerChannel;
    srate=SDATA.info.sampling_rate;

    for c=1:size(trigcodes,1) % For conditions

        % Segmentation
        [alpha_phase, isNotArtifact, timeVecITPC]=segmentContEEGdata(trigcodes{c} , time,...
            alpha_phase_whole_ts , triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        artrej=1;
        if artrej==1
            Alpha_SingleTrials{c}=alpha_phase(:,:,logical(isNotArtifact));
        end

        ITPCalpha_nTrials(c)=size(Alpha_SingleTrials{c},3);
        ITPCalpha_timevec{c}=timeVecITPC;

        % Save
        cd  'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Alpha Phase'
        save(savefilename, "Alpha_SingleTrials","ITPCalpha_nTrials","ITPCalpha_timevec",'-v7.3')

    end
end
end