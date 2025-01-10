%% Change across time window in instantaneous frequency

clear
clc
subj=[3];
elec=[25:30 62:64]; % occipital electrodes
srate=1024; % sampling rate in hertz

% Run
cd  'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Alpha Phase'
for s=1:length(subj)
    figure; t=tiledlayout('flow');   
    maintit=sprintf("Subject %i",subj(s));
    title(t,maintit)
    for t=1:2 % for the two target conditions (target and catch)
        if t==1
            % Load Target Trials
            loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Target.mat',subj(s));
            load(loadfilename1, "Alpha_SingleTrials","ITPCalpha_timevec")
        elseif t==2
            % Load Target Trials
            loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Catch.mat',subj(s));
            load(loadfilename1, "Alpha_SingleTrials","ITPCalpha_timevec")
        end

    % Prep (separate by condition, select electrodes and average)
    cond1_angles=Alpha_SingleTrials{1,1};
    cond2_angles=Alpha_SingleTrials{1,2};

    cond1_angles=squeeze(circ_mean(cond1_angles(:,elec,:),[],2));
    cond2_angles=squeeze(circ_mean(cond2_angles(:,elec,:),[],2));

    % Average across trials
    cond1_angles=circ_mean(cond1_angles,[],2);
    cond2_angles=circ_mean(cond2_angles,[],2);

    % Calculate Instantaneous Frequency at each time point
    instalpha_cond1=diff(unwrap(cond1_angles))/(2*pi*srate);
    instalpha_cond2=diff(unwrap(cond2_angles))/(2*pi*srate);
    
    timeVec=ITPCalpha_timevec{1,1};

    % Apply suprathreshold filter (based on mike cohen video)
    instaz_1=(instalpha_cond1-mean(instalpha_cond1))/std(instalpha_cond1); % zscore
    histogram(instaz_1,200)
    instaz_2=(instalpha_cond2-mean(instalpha_cond2))/std(instalpha_cond2);
    histogram(instaz_2,200)

    tofilter_1=find(abs(instaz_1)>2); % find supra-threshold data points
    tofilter_2=find(abs(instaz_2)>2);

    % Replace those values with the median of the surrounding values
    instalpha_filt1=instalpha_cond1;
    instalpha_filt2=instalpha_cond2;
    k=round(1000*50/srate); % kernelsize
    for i=1:length(tofilter_1) % condition 1
        indices=max(1,tofilter_1(i)-k):min(length(instalpha_filt1),tofilter_1(i)+k);
        instalpha_filt1(tofilter_1(i))=median(instalpha_filt1(indices));
    end

    for i=1:length(tofilter_2) % condition 2
        indices=max(1,tofilter_2(i)-k):min(length(instalpha_filt2),tofilter_2(i)+k);
        instalpha_filt2(tofilter_2(i))=median(instalpha_filt2(indices));
    end

    % Plot Phase
    nexttile; 
    figure;plot(timeVec(1:end-1),instalpha_filt2);
   figure; plot(timeVec(1:end-1), instalpha_cond2)

    if t==1
        title("Target Trials 800")
    else
        title("Target Trials 850")
    end
    nexttile; plot(timeVec(1:end-1),instalpha_cond2);
    %nexttile;circ_plot(cond2_angles,'hist',[],20,false,false,'linewidth',2,'color','r')
    if t==1
        title("Catch Trials 800")
    else
        title("Catch Trials 850")
    end

%     % Save angles
%     angles(s,t,1)={cond1_angles}; % (target type, condition, trials)
%     angles(s,t,2)={cond2_angles};
    end
end

