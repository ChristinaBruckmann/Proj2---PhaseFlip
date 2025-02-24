%% Alpha Phase Angle Analysis

%% Plot Alpha ITPC Time Course

clear
clc

% What to plot?
subj=[1:9];
gl_plots=1; % also plot average across all participants in subj?
timecourse=0; % plot also the angle time course for each participant?

% Params 
elec=[25:30 62:64]; % occipital electrodes
timepoint=[800];
exclude_t=30; % how many trials excluded at the begining, so that learning effects can be minimized? (this currently happens after art rejection, so rejected trials are not equal across cond and subj, fix this!)

% Preallocate
GL_cond1_target=0;
GL_cond2_target=0;
GL_cond1_catch=0;
GL_cond2_catch=0;

% Dir
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
            % Load Catch Trials
            loadfilename1=sprintf('EEG_PF_Pilot_Subj%i_AlphaPhaseSegmented_Catch.mat',subj(s));
            load(loadfilename1, "Alpha_SingleTrials","ITPCalpha_timevec")
        end

    % Prep (separate by condition, select electrodes and average)
    cond1_angles=Alpha_SingleTrials{1,1};
    cond2_angles=Alpha_SingleTrials{1,2};

    cond1_angles=squeeze(circ_mean(cond1_angles(:,elec,:),[],2));
    cond2_angles=squeeze(circ_mean(cond2_angles(:,elec,:),[],2));

    % Select time point
    timeVec=ITPCalpha_timevec{1,1};
    [~,  timeidx] = min(abs(timeVec - timepoint));
    cond1_angles=cond1_angles(timeidx,exclude_t:end);
    cond2_angles=cond2_angles(timeidx,exclude_t:end);


    % Plot Angles
    nexttile; circ_plot(cond1_angles,'pretty','ro',true,'linewidth',2,'color','r');
    if t==1
        title("Target Trials 800")
    else
        title("Catch Trials 800")
    end
    nexttile; circ_plot(cond2_angles,'pretty','ro',true,'linewidth',2,'color','r');
    %nexttile;circ_plot(cond2_angles,'hist',[],20,false,false,'linewidth',2,'color','r')
    if t==1
        title("Target Trials 850")
    else
        title("Catch Trials 850")
    end

    % Save data for GL plots
    % Mean
    cond1_mean=circ_mean(cond1_angles);
    cond2_mean=circ_mean(cond2_angles);
    % Center trials around 0 from cond 1
    cond1_cen=cond1_angles-cond1_mean;
    cond2_cen=cond2_angles-cond1_mean; % IMPORTANT THAT THIS IS 1 AND NOT 2, so we can compare later!!

    if t==1
        GL_cond1_target=[GL_cond1_target cond1_cen];
        GL_cond2_target=[GL_cond2_target cond2_cen];
    else
        GL_cond1_catch=[GL_cond1_catch cond1_cen];
        GL_cond2_catch=[GL_cond2_catch cond2_cen];
    end
    end
end

%% Group Level Plots
if gl_plots
   maintitle=sprintf('Angles %i ms after WS', timepoint);
   figure; t=tiledlayout('flow');   title(t,maintitle)
   % Plot target trials
    nexttile;circ_plot(GL_cond1_target,'hist',[],30,true,true,'linewidth',3,'color','r')
    title("Target Trials 800")
    nexttile;circ_plot(GL_cond2_target,'hist',[],30,true,true,'linewidth',3,'color','r')
    title("Target Trials 850")

   % Plot catch trials
    nexttile;circ_plot(GL_cond1_catch,'hist',[],30,true,true,'linewidth',3,'color','r')
    title("Catch Trials 800")
    nexttile;circ_plot(GL_cond2_catch,'hist',[],30,true,true,'linewidth',3,'color','r')
    title("Catch Trials 850")
end
%% Plot TIme Course of phase angles
if timecourse
    % Load Data
    elec=[25:30 62:64]; % occipital electrodes
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

            cond1_angles=squeeze(mean(cond1_angles(:,elec,:),2));
            cond2_angles=squeeze(mean(cond2_angles(:,elec,:),2));

            % Plot time course of angles for each trial and mean
            nexttile; % 800 condition
            timeVec=ITPCalpha_timevec{1,1};
            for trial=1:size(cond1_angles,2)
                plot(timeVec,cond1_angles(:,trial))
                hold on
            end
            cond1_mean=circ_mean(cond1_angles, [], 2); % Plot mean
            plot(timeVec,cond1_mean,'LineWidth',4,'Color','r')
            % Specify the x-value where you want the marker/flag
            x_marker = 800; % Example x value (adjust as needed)
            y_marker = interp1(timeVec, cond1_mean, x_marker); % Get corresponding y value
            text=sprintf('(%0.0f, %0.2f)', x_marker, y_marker);
            xline(800,'LineWidth',3,'Label',text,'LabelOrientation','horizontal')
            if t==1
                title('800ms Target')
            else
                title('800ms Catch')
            end

            nexttile; % 850 condition
            for trial=1:size(cond2_angles,2)
                plot(timeVec,cond2_angles(:,trial))
                hold on
            end
            cond2_mean=circ_mean(cond2_angles, [], 2); % Plot mean
            plot(timeVec,cond2_mean,'LineWidth',4,'Color','r')
            % Specify the x-value where you want the marker/flag
            x_marker = 800; % Example x value (adjust as needed)
            y_marker = interp1(timeVec, cond2_mean, x_marker); % Get corresponding y value
            text=sprintf('(%0.0f, %0.2f)', x_marker, y_marker);
            xline(800,'LineWidth',3,'Label',text,'LabelOrientation','horizontal')
            if t==1
                title('850ms Target')
            else
                title('850ms Catch')
            end
        end
    end
end