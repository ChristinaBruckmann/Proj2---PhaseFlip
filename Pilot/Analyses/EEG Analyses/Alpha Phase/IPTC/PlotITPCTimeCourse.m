%% Plot Alpha ITPC Time Course

clear
clc
subj=[3];
elec=[25:30 62:64]; % occipital electrodes
bl=1; % baseline to theoretical ITPC?

% Load Data
cd  'Y:\el-Christina\PhaseFlip\PF_Pilot\Results\Alpha Phase'
for s=1:length(subj)
    % Load Target Trials
    loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Target.mat',subj(s));
    load(loadfilename1, "Alpha_SingleTrials","ITPCalpha_nTrials","ITPCalpha_timevec")
    timeVec(s,1,:)=ITPCalpha_timevec{1, 1};
    alpha_data(s,1,1,:,:)=circ_r(Alpha_SingleTrials{1, 1},[], [], 3); % subject, target(1)-catch(2), 800(1)-850(2),time points, electrodes
    alpha_data(s,1,2,:,:)=circ_r(Alpha_SingleTrials{1, 2},[], [], 3);
    ntrials(s,1,:)=ITPCalpha_nTrials;

    % Load Catch Trials
    loadfilename2=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Catch.mat',subj(s));
    load(loadfilename2, "Alpha_SingleTrials","ITPCalpha_nTrials","ITPCalpha_timevec")
    timeVec(s,2,:)=ITPCalpha_timevec{1, 1};
    alpha_data(s,2,1,:,:)=circ_r(Alpha_SingleTrials{1, 1},[], [], 3); % subject, target(1)-catch(2), 800(1)-850(2),time points, electrodes
    alpha_data(s,2,2,:,:)=circ_r(Alpha_SingleTrials{1, 2},[], [], 3);
    ntrials(s,2,:)=ITPCalpha_nTrials;
end


% Average across selected electrodes
alpha_data=mean(alpha_data(:,:,:,:,elec),5);


% Baseline Correction
if bl
    % Calculate chance ITPC
    for s=1:length(subj)
        for t=1:2 % for target types (target - catch)
            for c=1:2 % for timing conditions (800-850)
                % Chance angles
                angles=2*pi*rand(ntrials(s,t,c),10000);
                angle_means(s,t,c)=mean(circ_r(angles));
                angle_std(s,t,c)=std(circ_r(angles));
            end
        end
    end

    % Baseline correct
    alpha_data=alpha_data-angle_means;
end

% Average across subjects
alpha_data=squeeze(mean(alpha_data,1)); % Across subjects - % result: target (1)/catch(2) - condition (800/850) - time points

% Plot in one plot
figure; 
colourvec1(1,:,:)=[0.4940 0.1840 0.5560 1;
    0.4660 0.6740 0.1880 1];
colourvec1(2,:,:)=[0.4940 0.1840 0.5560 0.5;
0.4660 0.6740 0.1880 0.5];
for  t=1:2 % two target types (target and catch)
    for c=1:2 % two timing conditions (800 and 850)
        plot(squeeze(timeVec(1,1,:)),squeeze(alpha_data(t,c,:)),"LineWidth",4,"Color",colourvec1(t,c,:))
        hold on;
    end
end

xline(0)
xline(800)
xline(850)
legend("Target 800","Target 850", "Catch 800", "Catch 850")


% Plot
colourvec2=[0.4940 0.1840 0.5560 1;
    0.4660 0.6740 0.1880 1];

figure;tiledlayout('flow');
for  t=1:2 % two target types (target and catch)
    nexttile;
    for c=1:2 % two timing conditions (800 and 850)
        plot(squeeze(timeVec(1,1,:)),squeeze(alpha_data(t,c,:)),"LineWidth",4,"Color",colourvec2(c,:))
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
