%% Experiment Code Project 2 - Phase Flip
% To do:

    % Minor Fixes:
        % Recalculate max-frames
        % Option to adjust gabor intensity after block 1?

    % Major Steps:
        % Write Instructions
        % Test code functionality for the whole experiment

    % Discuss with Assaf:
        % Constraints for Matrix Shuffle
        % Add block label to the experiment?
        % Record also JND Task EEG?

% Requires the custom functions  PTB_plotfig.m and navipage.m, unlock_continue.m and respfunction.m

% Clear
sca;
close all;
clear;
% clear mex;

% Initialize Structs
scr=[]; % Everything related to PTB Screen
stim=[]; % Stimulus Information
time=[]; % Timing Information

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj2 - PhaseFlip\Pilot\Experiment Code Pilot'
%% Input GUI
input_complete=0;
reload_data=0; % default: do not reload
scr.unlock_code=["123"]; % code entered by experimenter before code continues (at certain points)

while ~ input_complete
    % Which type of settings would you like?
    title = 'What would you like to do?';
    options = {'Run Participant', 'Trigger Check', 'Oscilloscope Test','Debugging'};

    run_type = centeredMenu(title, options{:});

    switch run_type
        case 0
            %             f=errordlg("Please select what you would like to do.");
            %             uiwait(f)
            error('Terminated by user.')
        case 1 % Running a participant
            floatwin=1; % Take over whole screen (0) or just part of it (1)
            SkipSync=1; % Skip Sync Tests (1 when testing on laptop)
            testing=0; % testing the code? reduces amount of trials to a minimum
            speedrun=1; % automatically chooses a response, no need to manually click anything (used for testing the code)
            gamification=1; % let participant collect points?
            scr.scope=0; % run scope test?
            reminders=1; % Show reminder intervals?

            % Experiment Sections
            JND_task=0; % JND Task?
            practice=0; % Practice?
            staircase=0; % Staircase?

            % Input Subject Numbert
            sub_num_compl=0;
            subresults.subj_num=[];
            while ~sub_num_compl
                valid_subnumber=0;
                while ~valid_subnumber
                    subresults.subj_num =  str2double(string(inputdlg('Subject Number:')));
                    if ~isempty(subresults.subj_num) % If number has been entered, check if it is a valid number
                        isWholePositiveInteger = isnumeric(subresults.subj_num) && isfinite(subresults.subj_num) && subresults.subj_num > 0 && mod(subresults.subj_num,1) == 0;
                        if isWholePositiveInteger
                            valid_subnumber=1;
                        else
                            f=errordlg("Subject number must be a whole positive integer.");
                            uiwait(f)
                        end
                    else
                         error('Terminated by user.')
                    end
                end

                % Check if a file for this participant already exists
                savefilename=sprintf("Pilot_PhaseFlip_Subj%i",subresults.subj_num);
                try
                    load(savefilename) % Try to load a file with this name, if it exists ask if this should be used
                    answer = questdlg('Subject already exists. Continue?','Existing Subject', 'Yes','No, change number.','No, change number.');
                        switch answer
                            case 'Yes'
                                % Start where the participant has left off?
                                restart_text=sprintf('Reloading file.\n\nLast completed block: %i \n\nContinue with block %i?',subresults.status.last_block,subresults.status.last_block+1);
                                answer = questdlg(restart_text,'Existing Subject', 'Yes','No - restart.','Yes');
                                if ~isempty(answer)
                                    switch answer
                                        case 'Yes'
                                            % Mark the reload
                                            reload_data=1;
                                            sub_num_compl=1;
                                        case 'No - restart.'
                                            old_subresults=subresults; % save previous results in file under new name
                                            subresults=[]; % restart with clean subresults
                                            subresults.old_subresults=old_subresults;
                                            save(savefilename,"subresults")
                                            f=msgbox(sprintf("Previous results saved.\n\nNew file created."));
                                            uiwait(f)
                                            sub_num_compl=1;
                                    end
                                else
                                    error('Terminated by user.')
                                end
                            case 'No, change number.'
                                continue % Re-enter subject number
                        end
                catch ME
                    answer = questdlg('New subject. Create a new file?','New Subject', 'Yes','Return','Return');
                    % Handle response
                    if ~ isempty(answer)
                        switch answer
                            case 'Yes'
                                save(savefilename,"subresults");
                                f=msgbox("New file created.");
                                uiwait(f)
                                sub_num_compl=1;
                            case 'Return'
                                continue % Re-enter subject number
                        end
                    else
                        error('Terminated by user.')
                    end
                end
            end

            % Input Experimenter Name
            subresults.experimenter = {};
            while true
                answer = inputdlg('Experimenter Name:');
                if isempty(answer)  % User closed the dialog without entering anything
                    error('Terminated by user.')
                elseif ~isempty(strtrim(answer{1}))  % Check if input is non-empty and not just spaces
                    subresults.experimenter = string(answer);  % Valid input, exit the loop
                    break;
                end
            end

            input_complete=1;
        case 2 % Trigger Check
            % Jumps directly to the real experiment to test a trial, this skips JND, Staircase and Practice
            floatwin=0;
            SkipSync=0;
            testing=1;
            speedrun=0;
            scr.scope=0;
            gamification=0;
           reminders=0; % Show reminder intervals?

            % Experiment Sections
            JND_task=0;
            staircase=0;
            practice=0; % Easy practice?

            savefilename=sprintf('trigger_test_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            input_complete=1;
        case 3 % Oscilloscope Test
            floatwin=0;
            SkipSync=0;
            testing=1;
            speedrun=0;
            scr.scope=1;
            gamification=0;
            reminders=0; % Show reminder intervals?

            input_complete=1;
        case 4 % Debugging
            % Choose which debugging options you want
            prompt = 'Select debugging options:';
            options = {'Reduced_TrialN', 'Speedrun', 'JNDTask','Staircase','Practice','SkipSync','Floating_Window','Gamification','ReminderIntervals'};

            % Display list dialog and get selected indices
            [selectedIndices, ~] = listdlg('PromptString', prompt, ...
                'ListString', options, ...
                'SelectionMode', 'multiple', ...
                'OKString', 'Select', ...
                'CancelString', 'Cancel');

            % Initialize struct to store variables
            varStruct = struct();

            % Assign values based on selection
            for i = 1:length(options)
                if ismember(i, selectedIndices)
                    varStruct.(options{i}) = 1;
                else
                    varStruct.(options{i}) = 0;
                end
            end

            % Debugging Options Chosen
            testing = varStruct.(options{1});
            speedrun = varStruct.(options{2});
            JND_task = varStruct.(options{3});
            staircase = varStruct.(options{4});
            practice = varStruct.(options{5});
            SkipSync=varStruct.(options{6});
            floatwin=varStruct.(options{7});
            gamification=varStruct.(options{8});
            reminders=varStruct.(options{9});

            savefilename=sprintf('test_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            scr.scope=0;
            input_complete=1;
            subresults.subj_num=9999;
    end
end

%% PTB Set-Up

% Setup PTB
PsychDefaultSetup(2);
Priority(1);
Screen('Preference', 'SkipSyncTests', SkipSync);
scr.number = max(Screen('Screens'));

% Screen Values
white = WhiteIndex(scr.number);
scr.black = BlackIndex(scr.number);
scr.grey = white / 2;
if ~scr.scope
    scr.background = scr.grey;
    scr.fontcolour = scr.black;
else
    scr.background = scr.black;
    scr.fontcolour = white;
end

% Get screen size of main display
[scr.screenXpixels, scr.screenYpixels] = Screen('WindowSize', scr.number);

% Window size for the floating window
winWidth = 1600;
winHeight = 1600;

% Calculate center position for floating window
xPos = (scr.screenXpixels - winWidth) / 2;
yPos = (scr.screenYpixels - winHeight) / 2;

% Open Screen (with or without floating window)
if ~floatwin
    [scr.win, scr.windowRect] = PsychImaging('OpenWindow', scr.number, scr.background, [], 32, 2, [], [], kPsychNeed32BPCFloat);
else
    [scr.win, scr.windowRect] = PsychImaging('OpenWindow', scr.number, scr.background, ...
        [xPos, yPos, xPos + winWidth, yPos + winHeight], [], [], [], [], [], kPsychGUIWindow);
end

scr.ifi = Screen('GetFlipInterval', scr.win);

% Define Screen Dimensions
[scr.xCenter, scr.yCenter] = RectCenter(scr.windowRect);
[scr.axisx] = scr.windowRect([1, 3]);
[scr.axisy] = scr.windowRect([2, 4]);

% % Trigger Set-Up
% triggerPortAddress=hex2dec('FFF8');
% triggerPort=io64;
% s=io64(triggerPort) % Do not suppress output
% pause(1) % Add a pause so you can inspect the output
% io64(triggerPort, triggerPortAddress, 0); % Sets the trigger to zero for now


%% Timing Parameters

time.ITI=[1.85:0.1:2.15]; % Inter-Trial Interval
time.initmask=[0.75:0.1:1.25];% Initial Mask pre-WS
time.preq=[0.4:0.05:0.6]; % Mask Duration post target and pre-question
time.cuedur=0.1; % Cue Duration
time.tardur=0.016; % Target Duration
time.practiceinterval=0.7; % how long is the interval during practice? (not converted into frames here, happens in trial function)

% Target times for each block type / condition
time.trialtimes1=[0.7, 0.75, 1.2, 0.5, NaN]; %[timereg timeirr timelong timeshort timecatch];
time.trialtimes2=[0.75, 0.7, 1.25, 0.55, NaN]; %[timereg timeirr timelong timeshort timecatch];

% Save the timing as seconds before converting into frames
subresults.timeinsec=time;

% Convert into frames
time.ITI=round(time.ITI/scr.ifi); % Inter-Trial Interval
time.initmask=round(time.initmask/scr.ifi); % Initial Mask pre-WS
time.preq=round(time.preq/scr.ifi); % Pre-Question Interval
time.cuedur=round(time.cuedur/scr.ifi); % Cue Duration
time.tardur=round(time.tardur/scr.ifi); % Target Duration

% Others
time.maxresptime=3; % After which response delay is a warning displayed?

%% Stimulus Parameters

% Noise Parameters
stim.rectSize=300;
stim.noiseMean= 50;

% Gabor Parameters
stim.maskintensity=0.6;
practicegabor=0.7; % practice gabor intensity for easy practice

if reload_data && subresults.status.last_block>0
    gaborpercent=subresults.status.gaborintensity;
else % load previous gabor intensity
    gaborpercent=1-stim.maskintensity; % Initialize like this, adjust with staircase later
end

% Cue parameters for scope test
stim.baseRect = [0 0 120 120];
stim.rectColorCue = 1;

sigma = 0.4;
numCycles = 21; % keep an odd number?
degrees=pi/180; % units
theta=45*degrees; % tilt

% Fixation Cross
fixCrossDimPix = 10;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
stim.fixCoords = [xCoords; yCoords];
stim.lineWidthPix = 2;

% Cue Parameters
csLum=0; % luminance cue signal
circleweak= 0.5; % how weak is the circle compared to mask (lower means less contrast)

%% Text Parameters
text.keyleft="l";
text.keyright="r";

text.instructions={'During this experiment you will be presented with a  visual stimulus whose orientation you are asked to judge. \n The visual stimulus is an oriented Gabor patch, called the target.\n\n\n\n Continue with the arrow.';
     'Before the Gabor patch, you will see a grey circle. \n This is called the cue and it indicates that the target is about to be presented. \n\n The target (Gabor patch) will be left or right oriented. \n However, note that a target will not always be presented. No question will be asked, if that is case.'
     'During the experiment you will also asked to remember an indicated interval.\n The interval refers to the difference in time between the presentation of the first cue and a second cue. \n Within this interval, no target is shown.';
     'Before the start of the experimental blocks, your threshold will be established. \n This consits in presenting the target at different intensity-levels in order to find the optimal one for you! \n'
     'Lastly, during the experiment you will get the opportunity to obtain points. \n You will be awarded +1 point for every correct answer. The more points, the better your performance! \n\n And a good performance should be rewarded, right? And so it shall, indeed! \n\n\n\nIf you have understood the task and instructions, you can confirm and continue with the arrow.'};
text.JND_instructions={'During this task you will be asked to compare two intervals. \n The interval refers to the difference in time between the presentation of two cues.\n The target is not shown during the interval. \n\n\n\n Continue with the arrow.';
     'The intervals shown to you will be numbered. \n If the first interval is longer, then you must press 1. \n If the second interval is longer, then you must press 2. ';
     'Please note that it is important you do not move during this task, as any movement will introduce noise in the data. \n\n\n\n If you have understood the task and instructions, you can confirm and continue with the arrow.'};
text.short_instructions={'The visual stimulus that is presented to you is an oriented Gabor patch. \n This is called the target. \n\n Before the Gabor patch, you will see a grey circle. \n This is called the cue and it indicates that the target is about to be presented.';
     'The target (Gabor patch) will be left or right oriented. \n Press the left key if the Gabor patch is left-oriented.\n Press the right key if the Gabor patch is right-oriented. \n\n Please note that target will not be shown always! \n When the target is not shown, you do not have to press the left or right key.'};
text.easypractice_instructions={'Easy Practice. \n\n This is the part where you get to practice discriminating the orientation of the target. \n Remember that if you focus on the upper half of the target, it will be easier to discriminate the orientation. \n Press the left key if the Gabor patch is left-oriented.\n Press the right key if the Gabor patch is right-oriented. \n\n\n\n If you are ready, continue with the arrow.'};
text.practice_instructions={'Now you get to do a longer practice. \n\n Note that the cue will be shown first. \n The target will be presented after the cue is shown. \n\n Remember that if you focus on the upper half of the target, it will be easier to discriminate the orientation. \n\n\n If you are ready, continue with the arrow.'};
text.staircase_instructions={'Staircase. \n\n In order for you to achieve optimal performance during the experiment, we will first need to establish what is known as threshold. \n\n\n\n Continue with the arrow.';
    'The target will be shown to you at different contrast levels. \n This is done in order calculate the best intensity-level for you. \n At times it may be slighlty more difficult to discriminate the target, and the orientation, but\n do not worry because this is by design. \n If you are unsure of the orientation, venture your best guess!';
    'Please report the perceived orientation by using the left or right key. \n Remember: \n the target will always be shown, \n it is easier to discriminate the orientation if you focus on the upper half of the target. \n\n\n\n If you are ready, continue with the arrow.'};
text.block1_instructions={'Block 1. \n\n Use the arrows to continue';
     'Press the left key if the target (Gabor patch) is left-oriented.\n Press the right key if the target is right-oriented. \nWhen the target is not shown, you do not have to press the left or right key'};
text.block2_instructions={'Block 2.  \n\n Use the arrows to continue';
    'Press the left key if the target (Gabor patch) is left-oriented.\n Press the right key if the target is right-oriented. \nWhen the target is not shown, you do not have to press the left or right key'};
%% Create Trial Matrix and Design
% 400 trials in total (bit less than 45 min at 5,5 seconds per trial)
% 50 trials per block
% 8 blocks
if testing % shorter version to test the code
    % Blocks, Trials etc.
    cond=[1 2]; % conditions (1=800ms, 2=850ms)
    nblocks=length(cond); % total blocks
    ntrials=5; % trials per block
    ntrialstot=nblocks*ntrials;

    % Split per block
    nregular=1;
    nirregular=1;
    nlong=1;
    nshort=1;
    ncatchtrials=1;
else % experiment settings
    % Blocks, Trials etc.
    if subresults.subj_num % counterbalance
        cond=[1 2 1 2 1 2 1 2]; % conditions (1=800ms, 2=850ms)
    else
        cond=[2 1 2 1 2 1 2 1]; % conditions (1=800ms, 2=850ms)
    end
    nblocks=length(cond); % total blocks
    ntrials=50; % trials per block
    ntrialstot=nblocks*ntrials;

    % Split per block
    nregular=25;
    nirregular=0;
    nlong=4;
    nshort=4;

    % Left overs are assigned to catch trials (throw error if none are left)
    if (ntrials-nregular-nirregular-nlong-nshort)>49
        error('No catch trials left. Double check trial matrix generation. Catch trials are the most important part of this project.')
    elseif (ntrials-nregular-nirregular-nlong-nshort)>45
        warning('Using less than 5 catch trials per block. This is not advised. Please reconsider.')
    else
        ncatchtrials=ntrials-nregular-nirregular-nlong-nshort;
    end
end

if ~reload_data % only create new matrix if there is no previous one
    % Create trial matrix
    trialmatrix=[]; % Preallocate matrix
    for bl=1:nblocks % Create all blocks
        if cond(bl)==1
            trialtimes=time.trialtimes1;
        else
            trialtimes=time.trialtimes2;
        end
        timecatch=NaN; % here NaN (since invisible), in trial function random timing assigned for functionality

        ntrialtypes=[nregular nirregular nlong nshort ncatchtrials];

        % Create array of 50 trials with the desired split
        trialmatrix_temp=[];
        for tt=1:length(ntrialtypes)
            temp=[tt trialtimes(tt) cond(bl) bl]; % trialtype trialtime condition nblock
            temp=repmat(temp,ntrialtypes(tt),1); % repeat as many times as needed for this trial type
            trialmatrix_temp=[trialmatrix_temp; temp]; % append to block matrix
        end

        % Shuffle (ADD CONTRAINTS!)
        trialmatrix_temp=trialmatrix_temp(randperm(size(trialmatrix_temp,1)),:);

        % Append to main matrix
        trialmatrix=[trialmatrix; trialmatrix_temp];
    end
    subresults.trialmatrix=trialmatrix;
end

% %% Trigger Vector
% stim.triggervector=table();
%
% stim.triggervector(1,:)={100,101,110,111,112,113,114,120,121,122,130,131,140,141,142,143,150,151,152};
% stim.triggervector(2,:)={200,201,210,211,212,213,214,220,221,222,230,231,240,241,242,243,250,251,252};
% stim.triggervector(3,:)={1};% Set all to one for scope test
%
% stim.triggervector.Properties.VariableNames={'BlockStart','MaskOnset','CueOnsetRegular','CueOnsetIrregular','CueOnsetCatch','CueOnsetShort','CueOnsetLong','TargetOnsetLeft',...
%     'TargetOnsetRight','TargetOnsetCatch','MaskOffset','QuestionOnset','ResponseLeftCorrect', 'ResponseLeftIncorrect','ResponseRightCorrect','ResponseRightIncorrect','Reminder Interval Mask Onset','Reminder Interval Cue 1','Reminder Interval Cue 2'};
%% Create Stimuli

maxframesnoise=max(time.ITI)+max(time.ITI)+max(time.preq); % How many frames are needed per trial
maxframescue=time.cuedur;
DrawFormattedText(scr.win, 'Loading. Please wait.', 'center', 'center', scr.fontcolour);
Screen('Flip', scr.win);

%  Noise Mask
for fridx=1:maxframesnoise*4 % make 4 times the amount of frames just to be sure and also for JND task more are needed. FIX THIS!!!!
    if scr.scope % If oscillosope test, make mask as dark as possible
        stim.noiseimg=0.01+stim.maskintensity*(rand(stim.rectSize)-0.5);
        stim.noiseimg = min(max(stim.noiseimg, 0), 1); % confirm values are between 0 and 1
    else
        stim.noiseimg=0.5+stim.maskintensity*(rand(stim.rectSize)-0.5);
    end
    stim.noisetex(fridx)=Screen('MakeTexture', scr.win, stim.noiseimg); % Convert into Textures
end

% Cue Stimulus
[x, y]=meshgrid(-1:2/(stim.rectSize-1):1,-1:2/(stim.rectSize-1):1);
for fridx=1:maxframescue*2
    % Cue embedded in noise mask (from NoisyCircleCode.m)
    gausscw = exp(-1.*(x.^2+y.^2)./sigma.^2); % gaussian envelope for cue/warning
    gausscw=gausscw>0.3; % noisy but well defined circle (make bigger/smaller by changing this value)
    if ~scr.scope
        randMaskC=(rand(stim.rectSize)-0.5)*2; % Mask
        rescaledCS=(2*csLum-1)*gausscw; % Cue

        % Blend outside
        cueStim=rescaledCS;
        cueStim(gausscw<0.3)=0.5+0.5*(stim.maskintensity*randMaskC(gausscw<0.3)+(1-stim.maskintensity)*rescaledCS(gausscw<0.3));

        % Blend inside
        cueStim(gausscw>0.3)=0.5+0.5*(circleweak*randMaskC(gausscw>0.3)+(1-circleweak)*rescaledCS(gausscw>0.3));

        % colored cue (- and it's pink just for you)
        cueCol = repmat(cueStim, [1, 1, 3]); % expand the 2D matrix to a 3D one which is basically the color matrix
        cueCol(:, :, 1) = cueCol(:, :, 1) + 0.5*gausscw; % this is the red
        cueCol(:, :, 2) = cueCol(:, :, 2) + 0.2*gausscw; % this is green
        cueCol(:, :, 3) = cueCol(:, :, 3) + 0.3*gausscw; % this is blue
        stim.colCueTex(fridx) = Screen('MakeTexture', scr.win, cueCol);

    else
        % % Make Cues and Warning Signal Textures for scope test
        cueStim=0.01+stim.maskintensity*(rand(stim.rectSize)-0.5); % dark noisy mask
        xfrom=(length(cueStim)-stim.baseRect(3))/2; % Define cue rectangle
        xto=xfrom+stim.baseRect(3);
        cueStim(xfrom:xto,xfrom:xto)=stim.rectColorCue; % White cue, no blending for optimal luminance contrast
        cueStim = min(max(cueStim, 0), 1); % confirm values are between 0 and 1
    end
    stim.cuetex(fridx)=Screen('MakeTexture', scr.win, cueStim);  % Convert into Textures
end

% Gabor Patches
kxright=(pi*numCycles)*cos(theta);
kyright=(pi*numCycles)*sin(theta);
kxleft=(pi*numCycles)*cos(-theta);
kyleft=(pi*numCycles)*sin(-theta);

imRright=cos(kxright*x+kyright*y); % Create both orientations
imRleft= cos(kxleft*x+kyleft*y);

gaussEnvt = exp(-1.*(x.^2+y.^2)./sigma.^2); % create gaussian envelope (target)
gaussEnvt(gaussEnvt<0.01)=0; % change close-to-zero peripheral elements to zero
stim.gaborRight=imRright.*gaussEnvt;
stim.gaborLeft=imRleft.*gaussEnvt;
% Embedding it in mask happens lower in the script after the intensity has been chosen

%% Save all parameters and other info
if ~scr.scope
    if reload_data && subresults.status.last_block>0 % reload previous data
        resulttable=subresults.data; % otherwise initialize
    else
        subresults.status.JND_done=0;
        subresults.status.practice_done=0;
        subresults.status.staircase_done=0;
        subresults.status.last_block=0;
        subresults.status.total_points=0;
    end

    subresults.screeninfo=scr;
    subresults.textinfo=text;
    save(savefilename, "subresults")

    startblock=subresults.status.last_block+1;
end
%% Run Experiment
 if ~scr.scope
     try

         ListenChar(2); % Suppress key presses in command window

         % Introduction
         if ~reload_data
             starttext='Welcome to the experiment. \n\n Press any button to start.';
             navipage(scr.win, text.instructions) % long instructions
         else
             starttext='Welcome back. \n\n  Press any button to continue.';
         end

         DrawFormattedText(scr.win,starttext, 'center', 'center', scr.fontcolour);
         Screen('Flip', scr.win);
         KbStrokeWait;

         %% Run JND Task
         if ~ subresults.status.JND_done
             while JND_task
%                  io64(triggerPort, triggerPortAddress, 254); % Start recording
%                  WaitSecs(0.1)
%                  io64(triggerPort, triggerPortAddress, 0);
                 DrawFormattedText(scr.win,'Loading.', 'center', 'center', scr.fontcolour);
                 Screen('Flip', scr.win);
                 pause(3)

                 [JND, JND_results]=JNDfunction(scr,time,text.JND_instructions,stim);
                 % Plot results and flip to screen
                 Screen('Flip', scr.win);
                 JND_fig = figure('Visible', 'off'); % make invisible figure
                 plot(1:length(JND_results.Comparison), JND_results.Comparison, '-o', 'LineWidth', 2); xlim([1 length(JND_results.Comparison)]); ylim([min(JND_results.Comparison)-1 max(JND_results.Comparison)+1]); title('JND Results'); xlabel('Trial Number'); ylabel('Comparison Duration'); yline(JND) % plot data
                 PTB_plotfig(JND_fig, scr.win, "JND_Figure", 0) % plot figure to PTB screen with this custom function
                 Screen('Flip', scr.win);
                 pause(3)
%                  io64(triggerPort, triggerPortAddress, 255); % Stop recording
%                  WaitSecs(0.1)
%                  io64(triggerPort, triggerPortAddress, 0);

                 unlock_continue(scr.win, scr.unlock_code) % Blocks screen until experimenter unlocks (to prevent subject from changing the slide)

                 % Confirm or repeat
                 adjusttext=sprintf('Comparison Value: %f \n\n Confirm (1) or Repeat (0)?',JND);
                 [response]=respfunction(scr.win,adjusttext,["1","0"]);
                 if double(response)==0 % if adjustment requested, ask for confirmation
                     confirmationtext='Are you sure? \n\n Do not repeat unless you are sure something went wrong the first time. \n\n\n\n Accept staircase value (1) or Adjust anyway (0)?';
                     [response]=respfunction(scr.win,confirmationtext,["1","0"]);
                     if  double(response)==0
                         % Save previous results and repeat JND task
                         subresults.JND_res_discarded=JND_results;
                         subresults.JND_discarded=JND;
                         save(savefilename, 'subresults')
                     else
                         %Save results
                         subresults.JND_results=JND_results;
                         subresults.JND=JND;
                         save(savefilename, 'subresults')
                         JND_task=0; % Exit JND Task
                     end
                 else
                     %Save results
                     subresults.JND_results=JND_results;
                     subresults.JND=JND;
                     subresults.status.JND_done=1;
                     save(savefilename, 'subresults')
                     JND_task=0; % Exit JND Task
                 end
             end
         end

         %% Run Practice
         if practice
             easypractice=1;
         end

         while practice && ~subresults.status.practice_done
             navipage(scr.win,text.easypractice_instructions) % Show instructions
             % Easy Practice
             while easypractice
                 [~, RespEval, ~, ~]=trialfunction(scr,time,text,stim,practicegabor,[6, time.practiceinterval],speedrun);

                 % Print evaluation
                 if RespEval==1
                     DrawFormattedText(scr.win,'Correct', 'center', 'center', scr.fontcolour);
                     Screen('Flip', scr.win);
                     pause(1)
                 else
                     DrawFormattedText(scr.win,'Incorrect', 'center', 'center', scr.fontcolour);
                     Screen('Flip', scr.win);
                     pause(1)
                 end

                 % Another one? Or show instructions again? Or adjust gabor visibility?
                 while true
                     response=respfunction(scr.win,'Repeat? (1) \n\n Show Instructions (2)? \n\n Adjust Visibility? (3) \n\n Continue? (4)',["1","2","3","4"]);
                     switch double(response)
                         case 1
                             break % get out of the true loop and run another easy practice
                         case 2
                             navipage(scr.win,text.short_instructions) % show instructions
                             continue % then show the menu again
                         case 3
                             practicegabor=adjustintensity(scr,practicegabor); % adjust intensity
                             continue %before returning to the menu
                         case 4
                             easypractice=0; % finish easy practice
                             break % break out of true loop
                     end
                 end
             end
             subresults.status.practice_done=1;

             % Normal Practice
               navipage(scr.win,text.practice_instructions) % Show instructions
              [~, RespEval, ~, ~]=trialfunction(scr,time,text,stim,gaborpercent,[6, time.practiceinterval],speedrun);

                 % Print evaluation
                 if RespEval==1
                     DrawFormattedText(scr.win,'Correct', 'center', 'center', scr.fontcolour);
                     Screen('Flip', scr.win);
                     pause(1)
                 else
                     DrawFormattedText(scr.win,'Incorrect', 'center', 'center', scr.fontcolour);
                     Screen('Flip', scr.win);
                     pause(1)
                 end

                 % Another one? Or show instructions again? Or adjust gabor visibility?
                 while true
                     response=respfunction(scr.win,'Repeat? (1) \n\n Show Instructions (2)? \n\n Return to easy practice? (3) \n\n Continue? (4)',["1","2","3","4"]);
                     switch double(response)
                         case 1
                             break % get out of the true loop and run another easy practice
                         case 2
                             navipage(scr.win,text.short_instructions) % show instructions
                             continue % then show the menu again
                         case 3
                             easypractice=1; % turn on easy practice
                             break %before returning to the beginning
                         case 4
                             practice=0; % finish easy practice
                             break % break out of true loop
                     end
                 end

         end

         %% Run Staircase
         if staircase && ~subresults.status.staircase_done
             navipage(scr.win,text.staircase_instructions) % Show instructions
             [threshres,gaborpercent]=staircasefun(3,scr,time,text,stim,gaborpercent,speedrun);
         end

         % Accept Staircase Output or Adjust Gabor Difficulty
         if ~reload_data
             adjusttext=sprintf('Post-staircase intensity: %.2f \n\n Confirm (2) or Adjust (3)?',gaborpercent);
         else
             adjusttext=sprintf('Last Intensity: %.2f \n\n Confirm (2) or Adjust (3)?',gaborpercent);
         end

         [response]=respfunction(scr.win,adjusttext,["2","3"]);
         if double(response)==3 % if adjustment requested, ask for confirmation
             confirmationtext='Are you sure? \n\nIt is reccommended to stick to the staircase value. \n\n\n\n Accept staircase value (2) or Adjust anyway (3)?';
             [response]=respfunction(scr.win,confirmationtext,["2","3"]);
             if double(response)==3
                 gaborpercent=adjustintensity(scr, gaborpercent);
             end
         end

         % Save and Display the setting that will be used
         subresults.status.gaborintensity=gaborpercent;
         save(savefilename,"subresults")
         intconfirmed=sprintf('Confirmed Intensity: %.2f \n\nPress any button to continue.',gaborpercent);
         DrawFormattedText(scr.win,intconfirmed, 'center', 'center', scr.fontcolour);
         Screen('Flip', scr.win);
         KbStrokeWait;

         %% Run Blocks
         if reload_data
             tottrialcount=subresults.status.last_block*ntrials+1; % start at first trial of next block
         else
             tottrialcount=1; % total trial counter

             subresults.status.total_points=0;
             navipage(scr.win,text.block1_instructions) % Show instructions for first block
         end

         for b=startblock:nblocks
             starttext=sprintf('Ready? \n\n Press any key to start block %i',b);
             DrawFormattedText(scr.win,starttext, 'center', 'center', scr.fontcolour);
             Screen('Flip', scr.win);
             KbStrokeWait;

             %              io64(triggerPort, triggerPortAddress, 254); % Start recording
             %              WaitSecs(0.1)
             %              io64(triggerPort, triggerPortAddress, 0);
             DrawFormattedText(scr.win,'Loading.', 'center', 'center', scr.fontcolour);
             Screen('Flip', scr.win);
             pause(3)

             % Determine current condition
             currcond=cond(b);

             % Initialize gamification block points
             block_points=0;

             %              % Block Onset Trigger
             %              io64(triggerPort, triggerPortAddress,stim.triggervector{cond(b),'BlockStart'});
             %              WaitSecs(0.1)
             %              io64(triggerPort, triggerPortAddress, 0);

             % Run Trials
             for t=1:ntrials

                 trialinfo=subresults.trialmatrix(tottrialcount,:); % Choose trialinfo from trialmatrix

                 [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,gaborpercent,trialinfo,speedrun); % Run trial

                 if Resp~=9% If trial was not interruped, save and proceed as normal
                     resulttable(tottrialcount,:)=table(currcond, b, t, trialinfo(1), trialinfo(2),RT, Resp, RespEval, warning, gaborpercent,stim.maskintensity, 'VariableNames',{'Condition','Block', ...
                         'Trial','Trial Type','Target Interval','Reaction Time', 'Orientation Reseponse', 'Correct/Incorrect', 'Late Warning', 'Gabor Strength','Mask Intensity'});
                     subresults.data=resulttable;
                     tottrialcount=tottrialcount+1; % update total trial counter
                     %subresults.status.last_trial=tottrialcount; % save status
                     save(savefilename,"subresults");

                     % Gamification Update
                     if gamification && RespEval==1 % If correct, update points
                         subresults.status.total_points=subresults.status.total_points+1;
                         block_points=block_points+1;
                         DrawFormattedText(scr.win,'Correct! \n\n +1 Point', 'center', 'center', scr.fontcolour);
                         Screen('Flip', scr.win);
                         pause(0.5)
                     elseif gamification && RespEval==0
                         DrawFormattedText(scr.win,'Incorrect.', 'center', 'center', scr.fontcolour);
                         Screen('Flip', scr.win);
                         pause(0.5)
                     elseif gamification
                         DrawFormattedText(scr.win,'No target.', 'center', 'center', scr.fontcolour);
                         Screen('Flip', scr.win);
                         pause(0.5)
                     end

                     if reminders % If reminders are activated
                         if ~mod(t,4) % Show reminder interval every nth trial
                             timing_reminder(scr, time, stim, trialinfo)
                         end
                     end

                 else % If trial was interrupted, repeat the trial and do not save output
                     DrawFormattedText(scr.win,'Press any button to continue the task.', 'center', 'center', scr.fontcolour);
                     Screen('Flip', scr.win);
                     KbStrokeWait;
                 end
             end

             % Update Status
             subresults.status.last_block=b; % save status
             save(savefilename,"subresults");

             % Calculate block performance
             blockperf=resulttable{resulttable{:,'Block'}==b,'Correct/Incorrect'}; % performance results for latest block
             blockaccuracy=mean(blockperf,'omitnan'); % mean accuracy
             totalperf=resulttable{:,'Correct/Incorrect'};
             totalaccuracy= mean(totalperf,'omitnan'); % total accuracy of all trials SO FAR
             accuracy=[blockaccuracy totalaccuracy];

             DrawFormattedText(scr.win,'Loading.', 'center', 'center', scr.fontcolour);
             Screen('Flip', scr.win);
             pause(3)
             %              io64(triggerPort, triggerPortAddress, 255); % Stop recording
             %              WaitSecs(0.1)
             %              io64(triggerPort, triggerPortAddress, 0);

             % End of block/experiment message
             if b==nblocks
                 scoretext=sprintf('End of block %i/%i  \n\n Points in this block: %i \n\n Total score: %i ',b,nblocks,block_points, subresults.status.total_points);
                 DrawFormattedText(scr.win,scoretext, 'center', 'center', scr.fontcolour);
                 Screen('Flip', scr.win);
                 pause(2.5)

                 DrawFormattedText(scr.win,'This is the end of the experiment. \n\n Please wait for the experimenter.', 'center', 'center', scr.fontcolour);
                 Screen('Flip', scr.win);
                 unlock_continue(scr.win, scr.unlock_code) % Blocks screen until experimenter unlocks (to prevent subject from changing the slide)

                 DrawFormattedText(scr.win,'Data has been saved. \n\n Press any button to close the screen.', 'center', 'center', scr.fontcolour);
                 Screen('Flip', scr.win);
                 KbStrokeWait;
             else
                 next_block=0; % continue with next block?
                 while ~next_block
                     scoretext=sprintf('End of block %i/%i  \n\n Points in this block: %i \n\n Total score: %i ',b,nblocks,block_points, subresults.status.total_points);
                     DrawFormattedText(scr.win,scoretext, 'center', 'center', scr.fontcolour);
                     Screen('Flip', scr.win);
                     pause(2.5)
                     blockendmessage=sprintf('Please take a break. \n\n \n\nPress the space bar to continue with the next block.');
                     DrawFormattedText(scr.win,blockendmessage, 'center', 'center', scr.fontcolour);
                     % Print performance
                     performancetext=sprintf('B%.2f24753 /n T%.2f45264', accuracy(1),accuracy(2)); % block and total
                     [~,~,heightacc,lengthacc]=Screen('TextBounds', scr.win, performancetext); % get bounds
                     DrawFormattedText(scr.win, performancetext, scr.axisx(1), scr.axisy(2)-heightacc*1.2,scr.fontcolour);
                     Screen('Flip', scr.win);
                     [~, keyNamethr, ~]=KbStrokeWait;
                     % Exit?
                     if KbName(keyNamethr)=="space"
                         next_block=1;
                     elseif KbName(keyNamethr)=="e"
                         DrawFormattedText(scr.win,'Continue(c)? Exit(e)?', 'center', 'center',scr.fontcolour);
                         Screen('Flip', scr.win);
                         [~, keyNamethr, ~]=KbStrokeWait;
                         if KbName(keyNamethr)=="e"
                             ListenChar(0); % Restore default input handling
                             sca
                             return
                         elseif KbName(keyNamethr)=="c"
                             continue
                         end
                     end
                 end
             end
         end
         ListenChar(0); % Restore default input handling
         sca;
         ShowCursor;
     catch
         ListenChar(0); % Restore default input handling
         sca;
         ShowCursor;
         psychrethrow(psychlasterror);
     end
 else
     %% Scope Test
          scopedone=0;

     while  ~scopedone
     DrawFormattedText(scr.win,'Oscilloscope Test. \n\n Press any button to select timing.', 'center', 'center', scr.fontcolour);
     Screen('Flip', scr.win);
     KbStrokeWait;

     % Select interval duration
     adjusting=1;
     scope_dur=0.8;

     while adjusting
         adjustingtext=sprintf('Which interval do you want to test? \n\n\n Interval: %.2f seconds \n\n\n \n\n\nUp-Down with arrows. Confirm with enter. Exit with e.',scope_dur);
         DrawFormattedText(scr.win,adjustingtext, 'center', 'center', scr.fontcolour);
         Screen('Flip', scr.win);

         [keyDown, ~, keyName] = KbCheck;
         if keyDown
             if strcmp(KbName(keyName), 'UpArrow')==1
                 scope_dur=scope_dur+0.01;
                 pause(0.05)
             elseif strcmp(KbName(keyName), 'DownArrow')==1
                 scope_dur=scope_dur-0.01;
                 pause(0.05)
             elseif strcmp(KbName(keyName), 'Return')==1 % confirm with enter
                 newvalue=scope_dur;
                 adjusting=0;
             elseif strcmp(KbName(keyName),'e')==1 % exit with e
                 ListenChar(0); % Restore default input handling
                 scopedone=1;
                 sca
                 return
             end

             if scope_dur<0.02 % cannot go below 0.02 seconds
                 scope_dur=0.02;
             elseif scope_dur>4 % or above 4 seconds
                 scope_dur=4;
             end
         end
     end

     % Run trial function
     trialfunction(scr,time,text,stim,gaborpercent,[99,scope_dur,3],0); % Run trial

     end
 end
%% Trial Function
function [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,gaborpercentin,trialinfoin,speed)
% trialinfoin is the input for the trial information(type of trial and target interval duration from the trial matrix)
% Choose Inter stimulus interval based on trialinfo input and convert to frames
% speedrun automatically chooses a response so no manual response is needed (for code testing)

time.ISI=round(trialinfoin(2)/scr.ifi); % convert current target interval to frames

if trialinfoin(1)==5 % if current trial is a catch trial, make gabor invisible and assign random timing (since presentation time does not matter (but has to be significantly longer than 800ms to not smear in analysis)
    gaborpercentin=0;
    time.ISI=round(1.2/scr.ifi);
end

% Create gabors with the chosen difficulty
for i=1:time.tardur
    if scr.scope
        noiseimg=0.01+0.1*(rand(stim.rectSize)-0.5); % Dark noise mask
        % make cue area of noise tex white instead of showing Gabor
        xfrom=(length(stim.noiseimg)-stim.baseRect(3))/2;
        xto=xfrom+stim.baseRect(3);
        noiseimg(xfrom:xto,xfrom:xto)=stim.rectColorCue;
        stim.tartex_l(i)=Screen('MakeTexture', scr.win, noiseimg);
        stim.tartex_r=stim.tartex_l;
    else
        % Left
        randMask=(rand(stim.rectSize)-0.5)*2;
        targetStim=0.5+0.5*(stim.maskintensity*randMask+stim.gaborLeft*gaborpercentin);
        stim.tartex_l(i)=Screen('MakeTexture', scr.win, targetStim);
        % Right
        randMask=(rand(stim.rectSize)-0.5)*2;
        targetStim=0.5+0.5*(stim.maskintensity*randMask+stim.gaborRight*gaborpercentin);
        stim.tartex_r(i)=Screen('MakeTexture', scr.win, targetStim);
    end
end

% Shuffle textures and initialize counters
cnoise=1;
ccue=1;
cuetex=stim.cuetex(randperm(length(stim.cuetex)));
noisetex=stim.noisetex(randperm(length(stim.noisetex)));
tartex_l=stim.tartex_l(randperm(length(stim.tartex_l)));
tartex_r=stim.tartex_r(randperm(length(stim.tartex_r)));

% Randomly choose from jittered timings
ITI=time.ITI(randi(length(time.ITI)));
preq=time.preq(randi(length(time.preq)));
initmask=time.initmask(randi(length(time.initmask)));

% Choose orientation
if scr.scope % random for scope test
    targetorient=randi(2)-1;
    tartex=cuetex;
else % in theory based on trial matrix, but not yet implemented
    targetorient=randi(2)-1;
    if targetorient==1
        tartex=tartex_l; % Choose left stimuli
    else
        tartex=tartex_r; % Choose right stimuli
    end
end

% % Choose correct triggers for condition, trial type and orientation trvec=[CueOnset,TargetOnset]
%     if scr.scopetest
%         trvec=1;
%     else
%         % Which condition
%         switch trialinfoin(3)
%             case 1
%                 % Which Trial Type
%                 switch trialinfoin(1)
%                     case 1
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetRegular'}];
%                     case 2
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetIrregular'}];
%                     case 3
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetLong'}];
%                     case 4
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetShort'}];
%                     case 5
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetCatch'}];
%                 end
%
%                 % Which target orientation
%                 if trialinfoin(1)==5 % if catch trial
%                     trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetCatch'}];
%                 else % if not catch trial assign target orientation
%                     switch targetorient
%                         case 1
%                             trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetLeft'}];
%                         case 2
%                             trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetRight'}];
%                     end
%                 end
%             case 2
%                 % Which Trial Type
%                 switch trialinfoin(1)
%                     case 1
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetRegular'}];
%                     case 2
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetIrregular'}];
%                     case 3
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetLong'}];
%                     case 4
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetShort'}];
%                     case 5
%                         trvec=[stim.triggervector{trialinfoin(3),'CueOnsetCatch'}];
%                 end
%
%                 % Which target orientation
%                 if trialinfoin(1)==5 % if catch trial
%                     trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetCatch'}];
%                 else % if not catch trial assign target orientation
%                     switch targetorient
%                         case 1
%                             trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetLeft'}];
%                         case 2
%                             trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetRight'}];
%                     end
%                 end
%         end
%     end

% Start Trial

% Inter Trial Interval
for fixidx=1:ITI
    Screen('Flip', scr.win);
%     if fixidx==1
%         io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'MaskOnset'});
%     elseif fixidx==2
%         io64(triggerPort, triggerPortAddress, 0);
%     end
end

% Mask Onset
for fixidx=1:initmask
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
%     if fixidx==1
%         io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'MaskOnset'});
%     elseif fixidx==2
%         io64(triggerPort, triggerPortAddress, 0);
%     end
    cnoise=cnoise+1;
end

%Cue Onset
for cueidx=1:time.cuedur
    Screen('DrawTexture', scr.win, cuetex(ccue), [], [], 0);
    Screen('Flip', scr.win);
%     if cueidx==1
%         io64(triggerPort, triggerPortAddress,trvec(1)); % cue trigger
%     elseif cueidx==2
%         io64(triggerPort, triggerPortAddress, 0);
%     end
    ccue=ccue+1;
end

% Inter Stimulus Interval
for noiseidx=1:time.ISI
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
    cnoise=cnoise+1;
end

% Target Presentation
for taridx=1:time.tardur
    Screen('DrawTexture', scr.win,  tartex(taridx), [], [], 0);
    Screen('Flip', scr.win);
%     if taridx==1 % target trigger
%         io64(scr.triggerPort, scr.triggerPortAddress,trvec(2));
%     elseif taridx==2
%         io64(scr.triggerPort, scr.triggerPortAddress, 0);
%     end
end

% Post Target Mask
for noiseidx=1:preq
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
%     io64(scr.triggerPort, scr.triggerPortAddress, 0);
    cnoise=cnoise+1;
%     if noiseidx==preq
%         io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'MaskOffset'}); % Mask offset trigger
%     end
end


pause(0.01)
% io64(scr.triggerPort, scr.triggerPortAddress,0);

% Request Orientation Reponse
startTime=GetSecs;
responded=0;
currtime=0;
warning=0;
if  (trialinfoin(1)~=5) && ~speed && ~scr.scope % if not catch trial, in speed run mode or scope test, ask question as normal, otherwise automatically choose answer and continue
    while ~responded
        DrawFormattedText(scr.win, 'Left or right?', 'center', 'center', scr.fontcolour);
%         io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'QuestionOnset'});
%         pause(0.1)
%         io64(triggerPort, triggerPortAddress,0);
        if (currtime-startTime) > time.maxresptime % if participant hasn't responded after max time, display warning
            DrawFormattedText(scr.win,'Please answer.', 'center', scr.yCenter+250, scr.fontcolour);
            warning=1; % was a time warning displayed during response?
        end
        Screen('Flip', scr.win);
        [keyDown, respTime, keyName] = KbCheck;
        if keyDown
            if strcmp(KbName(keyName), text.keyleft)==1 % If left key
                RT=respTime-startTime;
                Resp=1; % left
                responded=1;
                if targetorient==1 % left
                    RespEval=1;
%                     io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseLeftCorrect'});
%                     pause(0.1)
%                     io64(triggerPort, triggerPortAddress,0);
                else
                    RespEval=0;
%                     io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseLeftIncorrect'});
%                     pause(0.1)
%                     io64(triggerPort, triggerPortAddress,0);
                end
            elseif strcmp(KbName(keyName), text.keyright)==1 % if right key
                RT=respTime-startTime;
                responded=1;
                Resp=0; %Right
                if targetorient==0 %Right
%                     io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseRightCorrect'});
%                     pause(0.1)
%                     io64(triggerPort, triggerPortAddress,0);
                    RespEval=1; % Evaluation
                else
%                     io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseRightIncorrect'});
%                     pause(0.1)
%                     io64(triggerPort, triggerPortAddress,0);
                    RespEval=0; % Evaluation
                end
            elseif strcmp(KbName(keyName),'e')==1
                DrawFormattedText(scr.win,'Continue (c)? Exit (e)?', 'center', 'center', scr.fontcolour);
                Screen('Flip', scr.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'e')
                    ListenChar(0); % Restore default input handling
                    sca
                    error("Terminated by user.")
                elseif strcmp(KbName(keyNamethr),'c')
                    Resp=9; % if trial was paused in between, do not go back to record late answers but continue with next trial instead
                    RespEval=0;
                    RT=9;
                    DrawFormattedText(scr.win,'Ready? Press any button.', 'center', 'center', scr.fontcolour);
                    Screen('Flip', scr.win);
                    KbStrokeWait
                    break
                end
            end
        end
        currtime=GetSecs;
    end

    % If the participant received a time out warning, remind them to respond faster next time
    if warning
        DrawFormattedText(scr.win,'Please remember to respond within 3 seconds.', 'center', 'center', scr.fontcolour);
        Screen('Flip', scr.win);
        WaitSecs(1.5)
    end
elseif trialinfoin(1)==5 % If catch, NaN in all fields
    RT=NaN;
    warning=NaN;
    Resp=NaN;
    RespEval=NaN;
else % in speed run assign random RT and random response.
    RT=rand(1);
    warning=0;
    Resp=randi([0,1]);
    if targetorient==Resp
        RespEval=1;
    else
        RespEval=0;
    end
end
end
%% Adjust Intensity Function
% Insert current intensity, adjust and return new intensity.
% Increase/Decrease with arrows, confirm with enter
function [newintensity]=adjustintensity(scr, gaborcontrast)
adjusting=1;
currlevel=gaborcontrast;
    while adjusting
    currentleveltext=sprintf('Current luminance level: %.2f',currlevel);
    DrawFormattedText(scr.win,currentleveltext, 'center', 'center', scr.fontcolour);
    Screen('Flip', scr.win);

    [keyDown, ~, keyName] = KbCheck;
    if keyDown
        if strcmp(KbName(keyName), 'UpArrow')==1
            currlevel=currlevel+0.01;
            pause(0.05)
        elseif strcmp(KbName(keyName), 'DownArrow')==1
            currlevel=currlevel-0.01;
            pause(0.05)
        elseif strcmp(KbName(keyName), 'Return')==1 % confirm with enter
            newintensity=currlevel;
            adjusting=0;
        elseif strcmp(KbName(keyName),'e')==1 % exit with e
            ListenChar(0); % Restore default input handling
            sca
            return
        end

        if currlevel<0.01 % cannot go below 0.01
            currlevel=0.01;
        elseif currlevel>0.99 % or above 0.99
            currlevel=0.99;
        end
    end
    end
end
%% Staircase Function
% Runs 1 up 3 down staircase
% Calls Trial Function
function [thresholdres,newgaborpercent]=staircasefun(reversals,scr,time,text,stim,oldgaborpercent, speedrun)
    finishstair=0;
    threshrun=1;
    %initstep=0.05; % initial step size adjustment until 1st rev
    stepsize=0.05; % after first reversal adjust by this step size
    currentgabor=oldgaborpercent;
    while ~finishstair
        thresholdfound=0;
        streak=0; % initiate correct response streak
        trialcount=0; % initiate trial counter
        reverscount=0; % initiate reversal counter
        firstrev=0; % initialize first reversal
        lastchange='undef' ; % initialize last change direction var
        while ~thresholdfound
            [resp, respeval, ~, ~]=trialfunction(scr,time,text,stim,currentgabor,[0, 0.725],speedrun); % Run trial function (last input is trial type (not relevant here) and target interval duration (between cond 1 and 2 here)r

            if ~isnan(resp) % if trial was not interrupted, continue as normal
                trialcount=trialcount+1;

                % Save
                thresholdres(threshrun,trialcount)=currentgabor;

                % Update
                if respeval==1 % if correct, increase streak
                    streak=streak+1;
                    if streak==3 % if third one in a row correct
                        currentgabor=currentgabor-stepsize; % decrease by step size
                        if strcmp(lastchange,'inc')
                            reverscount=reverscount+1; % update reversal counter
                        end
                        lastchange='dec'; % register last change as decrease
                        if firstrev==0 % If this was the first reversal, count the trial on which this happened
                            firstrev=trialcount;
                        end
                        streak=0; % restart streak counter
                        %prevstreak=1; % register that this was a streak
                    end
                else % if incorrect, reset streak and make task easier
                    streak=0;
                    %prevstreak=0; % reset previous streak variable
                    currentgabor=currentgabor+stepsize;
                    if strcmp(lastchange,'dec')
                        reverscount=reverscount+1; % update reversal counter
                    end
                    lastchange='inc';
                    if firstrev==0 % If this was the first reversal, count the trial on which this happened
                        firstrev=trialcount;
                    end
                end

                % End after n reversals has been reached
                if reverscount==reversals
                    thresholdfound=1;
                    newgaborpercent=mean(thresholdres(threshrun,firstrev:end)); % Calculate new threshold based on average from first reversal
                end

            else
                % if trial was interrupted, do not evaluate anything, do not update trial counter and just continue with the next iteration
                continue
            end
        end

        % Staircase Done. Please wait for experimenter
        DrawFormattedText(scr.win,'Task finished. Please wait for the experimenter.', 'center', 'center', scr.fontcolour);
        Screen('Flip', scr.win);
        KbStrokeWait;
        KbStrokeWait;
        KbStrokeWait;

        % Plot
        thresh_fig=figure('Visible', 'off');
        plot(thresholdres(threshrun,:));yline(newgaborpercent);ylim([0 1]); title('Staircase Result');xlabel('Trial');ylabel('Gabor Contrast')
        PTB_plotfig(thresh_fig, scr.win, "Staircase_Figure", 0) % plot figure to PTB screen with this custom function
        unlock_continue(scr.win, scr.unlock_code) % Blocks screen until experimenter unlocks (to prevent subject from changing the slide)

        % Repeat?
        [response]=respfunction(scr.win,'Repeat threshold? (y-2/n-3)',["2","3"]);
        if double(response)==2 % if repeat
            threshrun=threshrun+1;
        else % end staircase and save result
            finishstair=1;
        end

    end
end

%% Just-noticable-difference Function
% One up, one down for now. Change if needed
% Make sure to jitter the scr.background noise so that the longer interval is not obvious by noise alone.
% Maybe remove backfground noise?

function [JND,JND_results]=JNDfunction(scr,time,instr,stim)

navipage(scr.win,instr) % Show instructions

% Parameters
comparison_length=2; % start value for comparison interval length in s
standard_length=0.7; % length of standrad interval (ideally the same as the interval shown in the main task)
init_adjustment_steps=0.3; % bigger steps in the beginning (first descent)
adjustment_steps=0.05; % size of adjustment each step in seconds
maxreversals=5;
avg_trials=5; % average across how many of the last trials to calculate JND?

% Initialize
jndfound=0;
prevstep='down'; % preallocate previous response for reversal counter
revers_count=0; %preallocate reversal counter
JND_results=table(); %preallocate results
currtrial=0;

% Shuffle textures
cuetex=stim.cuetex(randperm(length(stim.cuetex)));
noisetex=stim.noisetex(randperm(length(stim.noisetex)));

while ~jndfound

    currtrial=currtrial+1; %Update trial counter
    cnoise=1;
    ccue1=1;
    ccue2=1;

    % Shuffle which interval gets presented first
    interval_lengths=[round(standard_length/scr.ifi), round(comparison_length/scr.ifi)]; % combine and convert into frames
    interval_lengths=interval_lengths(randperm(length(interval_lengths)));

    % Show Intervals
    for intervals=1:length(interval_lengths)
        % Fixation
        for fixidx=1:time.ITI
            Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour); % Adds current interval number to screen
            Screen('Flip', scr.win);
            cnoise=cnoise+1;
        end

        % Interval
        for cueidx=1:time.cuedur
            Screen('DrawTexture', scr.win, cuetex(ccue1), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
            Screen('Flip', scr.win);
            ccue1=ccue1+1;
        end

        for fixidx=1:interval_lengths(intervals)
            Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
            Screen('Flip', scr.win);
            cnoise=cnoise+1;
        end

        for cueidx=1:time.cuedur
            Screen('DrawTexture', scr.win, cuetex(ccue2), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
            Screen('Flip', scr.win);
            ccue2=ccue2+1;
        end

        % Fixation
        for fixidx=1:time.ITI
            Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
            Screen('Flip', scr.win);
            cnoise=cnoise+1;
        end

        % Blank screen to mark end of interval
        for fixidx=1:time.ITI
            Screen('Flip', scr.win);
        end
    end

    % Query Response
     [response]=respfunction(scr.win,'Which interval was longer? \n\n 1 or 2?',["1","2"]);

    % Evaluate Response
    if  interval_lengths(1)==round(comparison_length/scr.ifi) && comparison_length>standard_length  && double(response)==1 % if (comparison was presented first & is longer) AND (the answer was 1) --> correct
        resp_eval=1;
    elseif interval_lengths(2)==round(comparison_length/scr.ifi) && comparison_length>standard_length && double(response)==2 % if (comparison was presented second & is longer) AND (the answer was 2) --> correct
        resp_eval=1;
    elseif interval_lengths(2)==round(comparison_length/scr.ifi) && comparison_length==standard_length % standard and comparison are equal
        resp_eval=0; % always incorrect when comparison and standard are equal
    else
        resp_eval=0; % incorrect if none of the above hold
    end

    % Adjust Interval Length
    if revers_count==0 && resp_eval==1% until first reversal, go in bigger steps
        adj_step=init_adjustment_steps;
    else
        adj_step=adjustment_steps;
    end

    if  resp_eval==1 % response is correct
        % Lower comparison
        comparison_length_new=comparison_length-adj_step; % lower comparison toward standard
        currstep='down'; % register which direction the current step was taken
        if comparison_length_new<standard_length % if the new comparison would be smaller than the standard (can happen depending on which step size is chosen, so that it 'jumps over' the equal stage), make equal to standard
            comparison_length_new=standard_length;
        end
    else % if response is incorrect, make comparison longer
        comparison_length_new=comparison_length+adj_step;
        currstep='up'; % register which direction the current step was taken
    end

    % Reversal?
    if ~strcmp(currstep,prevstep) % if current direction is not equal to previous direction
        revers_count=revers_count+1; % count as reversal
    end
    prevstep=currstep; % update step direction for next iteration

        % Register data
        JND_results(currtrial,:)={standard_length, comparison_length, response, resp_eval,currstep};
    % Update Comparison or JND Found?
    if revers_count==maxreversals % if max reversals has been reached, take the average across the last 10 reversals as the JND
        jndfound=1;
        JND=mean(JND_results{end-avg_trials:end,2}); % average across last trials to establish JND
    else % if not found yet, update comparison
        comparison_length=comparison_length_new;
    end
end
JND_results.Properties.VariableNames={'Standard','Comparison','Response','Correct','Adjustment'}; % add variable names to results table
end

%% Timing reminder Function
% Reminder Intervals at full visibility
function timing_reminder(scr, time, stim, trialinfoin)

    % Randomly choose from jittered timings
    ITI=time.ITI(randi(length(time.ITI)));
    preq=time.preq(randi(length(time.preq)));
    initmask=time.initmask(randi(length(time.initmask)));
    time.ISI=round(trialinfoin(2)/scr.ifi); % convert current target interval to frames

    % Shuffle textures and initialize counters
    cnoise=1;
    ccue=1;
    cuetex=stim.colCueTex(randperm(length(stim.colCueTex)));
    noisetex=stim.noisetex(randperm(length(stim.noisetex)));

    % Text
    DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
    Screen('Flip', scr.win);
    pause (1)

    % Start Trial

    % Inter Trial Interval
    for fixidx=1:ITI
        DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
        Screen('Flip', scr.win);
    end

    % Mask Onset
    for fixidx=1:initmask
        DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
        Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
        Screen('Flip', scr.win);
        %     if fixidx==1
        %         io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'Reminder Interval Mask Onset'});
        %     elseif fixidx==2
        %         io64(triggerPort, triggerPortAddress, 0);
        %     end
        cnoise=cnoise+1;
    end

    %Cue Onset
    for cueidx=1:time.cuedur
        DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
        Screen('DrawTexture', scr.win, cuetex(ccue), [], [], 0);
        Screen('Flip', scr.win);
        %     if cueidx==1
        %         io64(triggerPort, triggerPortAddress,stim.triggervector{trialinfoin(3),'Reminder Interval Cue 1'}); % cue trigger
        %     elseif cueidx==2
        %         io64(triggerPort, triggerPortAddress, 0);
        %     end
        ccue=ccue+1;
    end

    % Inter Stimulus Interval
    for noiseidx=1:time.ISI
        DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
        Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0 );
        Screen('Flip', scr.win);
        cnoise=cnoise+1;
    end

    % Target Presentation
    for taridx=1:time.tardur
        DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
        Screen('DrawTexture', scr.win,  cuetex(ccue), [], [], 0);
        %     if taridx==1 % target trigger
        %         io64(triggerPort, stim.triggervector{trialinfoin(3),'Reminder Interval Cue 2'});
        %     elseif taridx==2
        %         io64(triggerPort, triggerPortAddress, 0);
        %     end
        Screen('Flip', scr.win);
    end

    % Post Target Mask
    for noiseidx=1:preq
        DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
        Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
        Screen('Flip', scr.win);
        cnoise=cnoise+1;
    end

end
