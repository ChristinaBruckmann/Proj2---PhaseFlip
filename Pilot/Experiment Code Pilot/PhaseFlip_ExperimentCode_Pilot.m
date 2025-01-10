%% Experiment Code Project 2 - Phase Flip

% Requires the custom functions PTB_plotfig.m and navipage.m, unlock_continue.m and respfunction.m, as well as palamedes toolbox

% Clear
sca;
close all;
clear;
clc;
% clear mex;

% Initialize Structs
scr=[]; % Everything related to PTB Screen
stim=[]; % Stimulus Information
time=[]; % Timing Information

%cd 'C:\Experiments\Christina\PhaseFlip\Pilot\FP_Data'
%% Input GUI
input_complete=0;
reload_data=0; % default: do not reload
scr.unlock_code=["123"]; % code entered by experimenter before code continues (at certain points)

while ~ input_complete
    % Which type of settings would you like?
    GUItitle = 'What would you like to do?';
    options = {'Run Participant', 'Trigger Check', 'Oscilloscope Test','Debugging'};

    run_type = centeredMenu(GUItitle, options{:});

    switch run_type
        case 0
            %             f=errordlg("Please select what you would like to do.");
            %             uiwait(f)
            error('Terminated by user.')
        case 1 % Running a participant
            floatwin=0; % Take over whole screen (0) or just part of it (1)
            SkipSync=0; % Skip Sync Tests (1 when testing on laptop)
            stim.reducetrials=0; % testing the code? reduces amount of trials to a minimum
            speedrun=0; % automatically chooses a response, no need to manually click anything (used for testing the code)
            gamification=1; % let participant collect points?
            scr.scope=0; % run scope test?
            reminders=1; % Show reminder intervals?

            % Experiment Sections
            JND_task=0; % JND Task?
            practice=1; % Practice?

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
            stim.reducetrials=1;
            speedrun=0;
            scr.scope=0;
            gamification=0;
            reminders=0; % Show reminder intervals?

            % Experiment Sections
            JND_task=0;
            practice=0; % Easy practice?

            savefilename=sprintf('trigger_test_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            input_complete=1;
        case 3 % Oscilloscope Test
            floatwin=0;
            SkipSync=0;
            stim.reducetrials=1;
            speedrun=0;
            scr.scope=1;
            gamification=0;
            reminders=0; % Show reminder intervals?

            input_complete=1;
        case 4 % Debugging
            % Choose which debugging options you want
            prompt = 'Select debugging options:';
            options = {'Reduced_TrialN', 'Speedrun', 'JNDTask','Practice','SkipSync','Floating_Window','Gamification','ReminderIntervals'};

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
            stim.reducetrials = varStruct.(options{1});
            speedrun = varStruct.(options{2});
            JND_task = varStruct.(options{3});
            practice = varStruct.(options{4});
            SkipSync=varStruct.(options{5});
            floatwin=varStruct.(options{6});
            gamification=varStruct.(options{7});
            reminders=varStruct.(options{8});

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

% Trigger Set-Up
scr.triggerPortAddress=hex2dec('FFF8');
scr.triggerPort=io64;
s=io64(scr.triggerPort) % Do not suppress output
pause(1) % Add a pause so you can inspect the output
io64(scr.triggerPort, scr.triggerPortAddress, 0); % Sets the trigger to zero for now


%% Timing Parameters

time.III_JND=[0.5:0.1:1]; % Inter-Interval-Interval for JND
time.ITI=[1.85:0.1:2.15]; % Inter-Trial Interval
time.ITIreminder=[0.85:0.1:1.15]; % ITI for the reminder function (shorter to not lose time)
time.initmask=[0.75:0.1:1.25];% Initial Mask pre-WS
time.preq=[0.4:0.05:0.6]; % Mask Duration post target and pre-question
time.cuedur=0.1; % Cue Duration
time.tardur=0.016; % Target Duration
time.practiceinterval=0.7; % how long is the interval during practice? (not converted into frames here, happens in trial function)
time.JNDpostintmask=[0.5:0.1:1.5]; % how long is the mask after each JND interval (jittered to not give away information)

% Target times for each block type / condition
time.trialtimes1=[0.7, 0.75, 1.2, 0.5, NaN]; %[timereg timeirr timelong timeshort timecatch];
time.trialtimes2=[0.75, 0.7, 1.25, 0.55, NaN]; %[timereg timeirr timelong timeshort timecatch];

% Save the timing as seconds before converting into frames
subresults.timeinsec=time;

% Convert into frames
time.III_JND=round(time.III_JND/scr.ifi); % Inter-Trial Interval
time.ITI=round(time.ITI/scr.ifi); % Inter-Trial Interval
time.ITIreminder=round(time.ITIreminder/scr.ifi); % Inter-Trial Interval
time.initmask=round(time.initmask/scr.ifi); % Initial Mask pre-WS
time.preq=round(time.preq/scr.ifi); % Pre-Question Interval
time.cuedur=round(time.cuedur/scr.ifi); % Cue Duration
time.tardur=round(time.tardur/scr.ifi); % Target Duration
time.JNDpostintmask=round(time.JNDpostintmask/scr.ifi); % JND mask

% Others
time.maxresptime=3; % After which response delay is a warning displayed?

%% Stimulus Parameters

% Noise Parameters
stim.rectSize=300;
stim.noiseMean= 50;

% Gabor Parameters
stim.maskintensity=0.6;
practicegabor=0.7; % practice gabor intensity for easy practice

% if reload_data && subresults.status.last_block>0
%     gaborpercent=subresults.status.gaborintensity;
% else % load previous gabor intensity
gaborpercent=1-stim.maskintensity; % Initialize like this, adjust with staircase later
% end

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

% % Load Images
[JNDimg, JNDColorMap] = imread("JND_illustration.png");
[intimg, intColorMap] = imread("Trial_Illustration.png");

% Turn into textures
JNDimgtex = Screen('MakeTexture', scr.win, JNDimg);
intimgtex = Screen('MakeTexture', scr.win, intimg);

%% Text Parameters
text.keyleft="l";
text.keyright="r";

text.instructions={'Welcome to the experiment! \n\n In this experiment we will test the limits of your timing skills.\n\n\n\n You can navigate the instructions with the arrows.';
    'In the first part, we will challenge you to time as precisely as possible.'};

text.JND_instructions1={'Part 1\n\n In this task, you will see two time-intervals.\n\n Your task is to say which one was longer. \n\n\n\n Continue with the arrow.';
    'Each interval is indicated by two black circles. \n\nPress any button to see what this looks like.\n\nYou do not have to respond yet.';};
text.JND_instructions2= {'Let`s practice this! \n\n\n\n Press the arrow to start the practice.'};
text.JND_instructions3= {'Great! \n\n You are ready to start the task. \n\n The task will become progressively more difficult. \n\n If you are not sure which interval was longer, just guess. \n\n\n\n Press the arrow to start the timing task.'};
text.JND_instructions4={'Well done, you have finished the timing test! \n\n\n\n Please wait for the experimenter.'};

text.maintask_instructions1={'Part 2 \n\n Now, we will start the main part of the experiment. \n\n In the following, you will see only one time-interval, very similar to those you saw before. \n\n\n\n Continue with the arrow.'};
text.maintask_instructions2={'Instead of a second black circle, we will show you stripes. \n\nYou task is to tell us, whether the stripes were tilted to the left or to the right. \n\nUse the l and r buttons for this.\n\n\n\n Continue with the arrow.' ;
    'We will now show you what this will look like. \n\n\n\n Continue with arrow.'};
text.maintask_instructions3={'Let`s practice this! \n\nPay close attention, this is very quick! \n\n\n\n Press the arrow to start the practice.'};
text.maintask_instructions4={['Great! Now that the task is clear, there are some important things you should know before we start: \n\n 1) The time between the black circle and the stripes is always the same. \n Learn this interval to focus at the right moment!\n This will help you see the target.' ...
    '\n\nWe will occasionally remind you of this interval!\n\n\n\n Continue with arrow.'];
    'One more thing: \n\n During the experiment you have the opportunity to collect points! \nYou will be awarded +1 point for every correct answer. \n\n And a good performance should be rewarded, right? \n The best 10% will receive a small additional payment!\n\n\n\n Continue with arrow.';
    'Let`s start with some easy targets! \n\n Remember to use the repeating timing to focus at the right moment! \n\n\n\n Press the arrow to start the task.'};
task.maintask_harder={'Now we will make the target harder to see. \n\n Sometimes you might not be able to see it - no worries, just take a guess! \n\n Every now and there will be no stripes. In that case, you will not be asked about the direction. \n\n\n\n Press the arrow to start the task.'};

text.block2_instructions={'Part 3.  \n\n Use the arrows to continue'; 'This part is almost identical to the previous one. \n\n The only difference is, that we slightly changed the time interval between the circle and the stripes.\n\n Pay close attention to learn the new interval.'; 'We will start with a few easy trials again. \n\n\n\n Start the task with the arrow. '};

%% Create Trial Matrix and Design
% 300 trials in total (bit less than 45 min at 5,5 seconds per trial)
% 50 trials per block

if stim.reducetrials % shorter version to test the code
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
    if mod(subresults.subj_num,2) % counterbalance
        cond=[1 2]; % conditions (1=800ms, 2=850ms)
    else
        cond=[2 1]; % conditions (1=800ms, 2=850ms)
    end
    nblocks=length(cond); % total blocks
    ntrials=150; % trials per block
    ntrialstot=nblocks*ntrials;

    % Split per block
    nregular=80;
    nirregular=0;
    nlong=12;
    nshort=12;

    % Left overs are assigned to catch trials (throw error if none are left)
    if (ntrials-nregular-nirregular-nlong-nshort)>=ntrials
        error('No catch trials left. Double check trial matrix generation. Catch trials are the most important part of this project.')
    elseif (ntrials-nregular-nirregular-nlong-nshort)>(ntrials-5)
        warning('Using less than 5 catch trials per block. This is not advised. Please reconsider.')
    else
        ncatchtrials=ntrials-nregular-nirregular-nlong-nshort;
    end
end

if ~reload_data % only create new matrix if there is no previous one
    % Create trial matrix (final matrix columns: trial type, trial timing, condition of the block, nblock)
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

        % Shuffle with constraints
        cond_fulfilled=0;
        while ~cond_fulfilled
            trialmatrix_temp=trialmatrix_temp(randperm(size(trialmatrix_temp,1)),:);

            if ~stim.reducetrials % only makes sense if trial number is not reduced
                % Rules
                first_regular=trialmatrix_temp(1:10,1)==1; % first x (10) trials need to be regular (for training purposes)
                catchConsec=max(diff([0; find(trialmatrix_temp(:,1)~=5); (size(trialmatrix_temp,1)+1)])-1); % maximum catch trials in a row

                if (mean(first_regular)==1) && (catchConsec<4) % if rules are fulfilled (first trials are all regular and no more than 3 consequtive catch trials in a row)
                    cond_fulfilled=1;
                end
            else
                cond_fulfilled=1;
            end
        end
        trialmatrix=[trialmatrix; trialmatrix_temp]; % concatenate trialmatrices for blocks
    end
    subresults.trialmatrix=trialmatrix;
end


%% Trigger Vectors

% Trigger Vector JND Task
stim.triggervectorJND=table();

stim.triggervectorJND(1,:)={007,008,009,010,011,012,030,031,032,033,034}; %Interval 1
stim.triggervectorJND(2,:)={007,008,009,020,021,022,030,031,032,033,034}; % Interval 2
stim.triggervectorJND.Properties.VariableNames={'JND Start','JND Stop','JND Trial Onset','JND Mask Onset','JND Cue 1','JND Cue 2', 'JND Question','JND Response 1 Correct','JND Response 1 Incorrect','JND Response 2 Correct', 'JND Response 2 Incorrect'};

% Trigger Vector Main Task
stim.triggervector=table();

stim.triggervector(1,:)={100,101,110,111,112,113,114,120,121,122,130,131,140,141,142,143,150,151,152};
stim.triggervector(2,:)={200,201,210,211,212,213,214,220,221,222,230,231,240,241,242,243,250,251,252};
stim.triggervector(3,:)={1};% Set all to one for scope test

stim.triggervector.Properties.VariableNames={'BlockStart','MaskOnset','CueOnsetRegular','CueOnsetIrregular','CueOnsetCatch','CueOnsetShort','CueOnsetLong','TargetOnsetLeft',...
    'TargetOnsetRight','TargetOnsetCatch','MaskOffset','QuestionOnset','ResponseLeftCorrect', 'ResponseLeftIncorrect','ResponseRightCorrect','ResponseRightIncorrect','Reminder Interval Mask Onset','Reminder Interval Cue 1','Reminder Interval Cue 2'};
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

% Create Gabor Example and text locations for subtitles
exgableft=[scr.axisx(2)*0.25, scr.axisy(2)*0.75]; % loc of gabor
exgabright=[scr.axisx(2)*0.75, scr.axisy(2)*0.75]; % loc of gabor
gaborRect = [0 0 180 180];
gaborRectleft = CenterRectOnPointd(gaborRect, exgableft(1),exgableft(2));
gaborRectright = CenterRectOnPointd(gaborRect, exgabright(1),exgabright(2));

if ~scr.scope
    exampletargetr = 0.5+0.5*(1*stim.gaborRight);
    exampletargetl = 0.5+0.5*(1*stim.gaborLeft);
else
    exampletargetr = 0.01+0.01*(1*stim.gaborRight);
    exampletargetl = 0.01+0.01*(1*stim.gaborLef);
end

exampletargetr = Screen('MakeTexture', scr.win, exampletargetr);
exampletargetl = Screen('MakeTexture', scr.win, exampletargetl);

[~,~,~,length1]=Screen('TextBounds', scr.win, 'Left'); % get bounds
[~,~,~,length2]=Screen('TextBounds', scr.win, 'Right'); % get bounds
textleft=[exgableft(1) round(exgableft(2)+gaborRect(4)*0.6)]; % shift on y axis in relation to gabor
textright=[exgabright(1) round(exgabright(2)+gaborRect(4)*0.6)]; % shift on y axis in relation to gabor
textleft(1)= textleft(1)-round(length1*0.5); % Centre text on location on x axis
textright(1)= textright(1)-round(length2*0.5); % Centre text on location on x axis
%% Save all parameters and other info
if ~scr.scope
    if reload_data && subresults.status.last_block>0 % reload previous maint task data if blocks have been done
        resulttable=subresults.data; % otherwise initialize
    elseif ~reload_data
        subresults.status.JND_done=0;
        subresults.status.practice_done=0;
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
            navipage(scr.win, text.instructions) % long instructions
        else
            starttext='Welcome back. \n\n  Press any button to continue.';
            DrawFormattedText(scr.win,starttext, 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            KbStrokeWait;
        end

        %% Run JND Task
        if ~ subresults.status.JND_done
            while JND_task

                % JND Introduction and Practice
                navipage(scr.win,text.JND_instructions1) % Show instructions part 1
                Screen('DrawTexture', scr.win, JNDimgtex);% show illustration
                Screen('Flip', scr.win);
                KbStrokeWait;
                navipage(scr.win,text.JND_instructions2) % Show instructions part 2

                % Practice JND (not yet implemented)
                practice_JND=1;
                while practice_JND

                    [ ~,resp_eval,~]=JNDfunction(scr,time,stim,1);
                    % Print evaluation
                    if resp_eval==1
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
                        response=respfunction(scr.win,'Repeat? (1) \n\n Show Instructions (2)? \n\n Continue? (3)',["1","2","3"]);
                        switch double(response)
                            case 1
                                break % get out of the true loop and run another easy practice
                            case 2
                                navipage(scr.win,text.JND_instructions1) % Show instructions part 1
                                Screen('DrawTexture', scr.win, JNDimgtex);% show illustration
                                Screen('Flip', scr.win);
                                KbStrokeWait;
                                navipage(scr.win,text.JND_instructions2) % Show instructions part 2
                                continue % then show the menu again
                            case 3
                                practice_JND=0; % finish easy practice
                                break % break out of true loop
                        end
                    end
                end

                % Main JND Task

                navipage(scr.win,text.JND_instructions3) % Show instructions part 3 
                io64(scr.triggerPort, scr.triggerPortAddress, 254); % Start recording
                WaitSecs(0.1)
                io64(scr.triggerPort, scr.triggerPortAddress, 0);

                DrawFormattedText(scr.win,'Loading.', 'center', 'center', scr.fontcolour);
                Screen('Flip', scr.win);
                pause(3)

                [JND_results,~,JND_UD]=JNDfunction(scr,time,stim,0);

                io64(scr.triggerPort, scr.triggerPortAddress, 255); % Stop recording
                WaitSecs(0.1)
                io64(scr.triggerPort, scr.triggerPortAddress, 0);

                % Save results
                cd 'C:\Experiments\Christina\PhaseFlip\Pilot\FP_Data'
                subresults.JND_results=JND_results;
                subresults.JND_UD=JND_UD;
                save(savefilename, 'subresults')

                % Calculate JND
                JND=mean(JND_UD.x(end-10:end));

                % Display end of JND message (wait for experimenter)
                navipage(scr.win,text.JND_instructions4)

                % Plot results and flip to screen
                t = 1:length(JND_UD.x); 
                JND_fig = figure('Visible', 'off'); % make invisible figure
                plot(t,JND_UD.x,'k'); hold on; plot(t(JND_UD.response == 1),JND_UD.x(JND_UD.response == 1),'ko', 'MarkerFaceColor','k'); plot(t(JND_UD.response == 0),JND_UD.x(JND_UD.response == 0),'ko', 'MarkerFaceColor','w'); set(gca,'FontSize',16); axis([0 max(t)+1 min(JND_UD.x)-(max(JND_UD.x)-min(JND_UD.x))/10 max(JND_UD.x)+(max(JND_UD.x)-min(JND_UD.x))/10]);
                JND_label=sprintf("JND=%.2f",JND); yline(JND); title(JND_label); xlabel('Trial'); ylabel('Stimulus Intensity'); 
                PTB_plotfig(JND_fig, scr.win, "JND_Figure", 0) % plot figure to PTB screen with this custom function
                Screen('Flip', scr.win);
                pause(3)


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
                        subresults.JND_UD_discarded=JND_UD;
                        save(savefilename, 'subresults')
                    else
                        %Save results
                        subresults.JND_results=JND_results;
                        subresults.JND_UD=JND_UD;
                        save(savefilename, 'subresults')
                        JND_task=0; % Exit JND Task
                    end
                else
                    %Save results
                    subresults.JND_results=JND_results;
                    subresults.JND_UD=JND_UD;
                    subresults.status.JND_done=1;
                    save(savefilename, 'subresults')
                    JND_task=0; % Exit JND Task
                end
            end
            subresults.status.JND_done=1;
        end

        %% Run Easy Practice
        if practice
            easypractice=1;
        end

        while practice && ~subresults.status.practice_done

            navipage(scr.win,text.maintask_instructions1) % Show instructions
            Screen('DrawTexture', scr.win, exampletargetr, [], gaborRectright);
            Screen('DrawTexture', scr.win, exampletargetl, [], gaborRectleft);
            DrawFormattedText(scr.win, 'Left', textleft(1),textleft(2),scr.fontcolour);
            DrawFormattedText(scr.win, 'Right', textright(1),textright(2),scr.fontcolour);
            navipage(scr.win,text.maintask_instructions2)
            Screen('DrawTexture', scr.win, intimgtex);% show illustration
            Screen('Flip', scr.win);
            KbStrokeWait;
            navipage(scr.win,text.maintask_instructions3)

            % Easy Practice
            while easypractice
                [~, RespEval, ~, ~]=trialfunction(scr,time,text,stim,practicegabor,[1, time.practiceinterval,1],speedrun);

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
                            navipage(scr.win,text.maintask_instructions1) % show instructions
                            navipage(scr.win,text.maintask_instructions4) % show instructions
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
        end

        %% Run Blocks

        % Reload data or initiate
        if reload_data
            tottrialcount=subresults.status.last_block*ntrials+1; % start at first trial of next block
        else
            tottrialcount=1; % total trial counter
            subresults.status.total_points=0;
            navipage(scr.win,text.maintask_instructions4) % Show instructions for first block
        end

        % Run blocks
        for b=startblock:nblocks

            % Preparation and Initialization
            currcond=cond(b); % Determine current condition
            block_points=0; % Initialize gamification block points

            % Initialize Palamedes Staircase Procedure
            up = 1;                     % inceasypracticerease after 1 wrong
            down = 3;                   % decrease after 3 consecutive right
            StepSizeDown = 0.05; % for the first trials, afterwards this gets reduced once we get closer to the threshold
            StepSizeUp = 0.05;
            stopcriterion = 'trials';
            stoprule = ntrials+1000; % needs to be defined, but I don't want it to stop by itself so I set it to higher than the trial number
            startvalue = gaborpercent;           %intensity on first trial
            xMax=1; % maximum gabor intensity
            xMin=0; % minimum gabor intensity

            UD = PAL_AMUD_setupUD('up',up,'down',down,'StepSizeDown',StepSizeDown,'StepSizeUp', ...
                StepSizeUp,'stopcriterion',stopcriterion,'stoprule',stoprule, ...
                'startvalue',startvalue, 'xMax',xMax,'xMin',xMin);

            % Show Instructions
            starttext=sprintf('Ready? \n\n Press any key to start block %i',b);
            DrawFormattedText(scr.win,starttext, 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            KbStrokeWait;

            io64(scr.triggerPort, scr.triggerPortAddress, 254); % Start recording
            WaitSecs(0.1)
            io64(scr.triggerPort, scr.triggerPortAddress, 0);
            DrawFormattedText(scr.win,'Loading.', 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            pause(3)

            % Block Onset Trigger
            io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{cond(b),'BlockStart'});
            WaitSecs(0.1)
            io64(scr.triggerPort, scr.triggerPortAddress, 0);

            % Run 20 Trials at full visbility
            for vistrials=1:20
                trialinfo=subresults.trialmatrix(tottrialcount,:); % Choose trialinfo from trialmatrix
                [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,UD.xCurrent,trialinfo,speedrun); % Run trial

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
            end

            navipage(scr.win, task.maintask_harder) % now we will make the tast harder - instructions

            % Run Real Trials
            for t=1:ntrials

                trialinfo=subresults.trialmatrix(tottrialcount,:); % Choose trialinfo from trialmatrix

                [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,UD.xCurrent,trialinfo,speedrun); % Run trial

                % Evaluate response and update everything

                if Resp~=9% If trial was not interruped, save and proceed as normal
                    resulttable(tottrialcount,:)=table(currcond, b, t, trialinfo(1), trialinfo(2),RT, Resp, RespEval, warning, UD.xCurrent,stim.maskintensity, 'VariableNames',{'Condition','Block', ...
                        'Trial','Trial Type','Target Interval','Reaction Time', 'Orientation Reseponse', 'Correct/Incorrect', 'Late Warning', 'Gabor Strength','Mask Intensity'}); % save trial in result table
                    subresults.data=resulttable; % save result table
                    tottrialcount=tottrialcount+1; % update total trial counter
                    save(savefilename,"subresults"); % save subresults

                    % Update Gabor Intensity (only if the trial was a regular one)
                    if trialinfo(1)==1
                        UD = PAL_AMUD_updateUD(UD, RespEval); % update UD structure
                    end

                    % Update step size after third reversal (to narrow the steps, assuming the participant has now reached a value around their threshold)
                    if ismember(3,UD.reversal)
                        UD = PAL_AMUD_setupUD(UD,'up',1,'down',1,'StepSizeDown',0.01,'StepSizeUp', ...
                            0.03);
                    end

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

                % Take a break every 50 trials (two breaks per 150 trials) but not at after the 150 trial
                if (~mod(t,50)) && (t~=ntrials)
                    io64(scr.triggerPort, scr.triggerPortAddress, 255); % Stop recording
                    WaitSecs(0.1)
                    io64(scr.triggerPort, scr.triggerPortAddress, 0);
                    DrawFormattedText(scr.win,'Please take a short break \n\n Press any button to continue', 'center', 'center', scr.fontcolour);
                    Screen('Flip', scr.win);
                    KbStrokeWait;
                    io64(scr.triggerPort, scr.triggerPortAddress, 254); % Start recording
                    WaitSecs(0.1)
                    io64(scr.triggerPort, scr.triggerPortAddress, 0);
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
            io64(scr.triggerPort, scr.triggerPortAddress, 255); % Stop recording
            WaitSecs(0.1)
            io64(scr.triggerPort, scr.triggerPortAddress, 0);

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
                    scoretext=sprintf('End of part %i/%i  \n\n Points in this part: %i \n\n Total score: %i ',b+1,nblocks+1,block_points, subresults.status.total_points);
                    DrawFormattedText(scr.win,scoretext, 'center', 'center', scr.fontcolour);
                    Screen('Flip', scr.win);
                    pause(3)
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
                        navipage(scr.win,text.block2_instructions)
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

% Choose correct triggers for condition, trial type and orientation trvec=[CueOnset,TargetOnset]
    if scr.scope
        trvec=1;
    else
        % Which condition
        switch trialinfoin(3)
            case 1
                % Which Trial Type
                switch trialinfoin(1)
                    case 1
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetRegular'}];
                    case 2
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetIrregular'}];
                    case 3
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetLong'}];
                    case 4
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetShort'}];
                    case 5
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetCatch'}];
                end

                % Which target orientation
                if trialinfoin(1)==5 % if catch trial
                    trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetCatch'}];
                else % if not catch trial assign target orientation
                    switch targetorient
                        case 1
                            trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetLeft'}];
                        case 0
                            trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetRight'}];
                    end
                end
            case 2
                % Which Trial Type
                switch trialinfoin(1)
                    case 1
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetRegular'}];
                    case 2
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetIrregular'}];
                    case 3
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetLong'}];
                    case 4
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetShort'}];
                    case 5
                        trvec=[stim.triggervector{trialinfoin(3),'CueOnsetCatch'}];
                end

                % Which target orientation
                if trialinfoin(1)==5 % if catch trial
                    trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetCatch'}];
                else % if not catch trial assign target orientation
                    switch targetorient
                        case 1
                            trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetLeft'}];
                        case 0
                            trvec=[trvec, stim.triggervector{trialinfoin(3),'TargetOnsetRight'}];
                    end
                end
        end
    end

% Start Trial

% Inter Trial Interval
for fixidx=1:ITI
    Screen('Flip', scr.win);
    if fixidx==1
        io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'MaskOnset'});
    elseif fixidx==2
        io64(scr.triggerPort, scr.triggerPortAddress, 0);
    end
end

% Mask Onset
for fixidx=1:initmask
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
    if fixidx==1
        io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'MaskOnset'});
    elseif fixidx==2
        io64(scr.triggerPort, scr.triggerPortAddress, 0);
    end
    cnoise=cnoise+1;
end

%Cue Onset
for cueidx=1:time.cuedur
    Screen('DrawTexture', scr.win, cuetex(ccue), [], [], 0);
    Screen('Flip', scr.win);
    if cueidx==1
        io64(scr.triggerPort, scr.triggerPortAddress,trvec(1)); % cue trigger
    elseif cueidx==2
        io64(scr.triggerPort, scr.triggerPortAddress, 0);
    end
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
    if taridx==1 % target trigger
        io64(scr.triggerPort, scr.triggerPortAddress,trvec(2));
    elseif taridx==2
        io64(scr.triggerPort, scr.triggerPortAddress, 0);
    end
end

% Post Target Mask
for noiseidx=1:preq
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
    cnoise=cnoise+1;
    if noiseidx==preq
        io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'MaskOffset'}); % Mask offset trigger
    elseif noiseidx==2
        io64(scr.triggerPort, scr.triggerPortAddress,0);
    end
end

Screen('Flip', scr.win);
pause(0.1)

% Request Orientation Reponse
startTime=GetSecs;
responded=0;
currtime=0;
warning=0;
if  (trialinfoin(1)~=5) && ~speed && ~scr.scope % if not catch trial, in speed run mode or scope test, ask question as normal, otherwise automatically choose answer and continue
    io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'QuestionOnset'});
    pause(0.1)
    io64(scr.triggerPort, scr.triggerPortAddress,0);
    while ~responded
        DrawFormattedText(scr.win, 'Left or right?', 'center', 'center', scr.fontcolour);
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
                    io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseLeftCorrect'});
                    pause(0.1)
                    io64(scr.triggerPort, scr.triggerPortAddress,0);
                else
                    RespEval=0;
                    io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseLeftIncorrect'});
                    pause(0.1)
                    io64(scr.triggerPort, scr.triggerPortAddress,0);
                end
            elseif strcmp(KbName(keyName), text.keyright)==1 % if right key
                RT=respTime-startTime;
                responded=1;
                Resp=0; %Right
                if targetorient==0 %Right
                    io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseRightCorrect'});
                    pause(0.1)
                    io64(scr.triggerPort, scr.triggerPortAddress,0);
                    RespEval=1; % Evaluation
                else
                    io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'ResponseRightIncorrect'});
                    pause(0.1)
                    io64(scr.triggerPort, scr.triggerPortAddress,0);
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
else % in speed run assign random RT and a probability response
    noise=-0.02 + (0.02 - (-0.02)) * rand; % noise of theoretical responder
    theothreshold = 0.5+noise; % threshold of theoretical responder
    RT=rand(1);
    warning=0;
    if theothreshold<=gaborpercentin % if gabor percent is above theoretical threshold count as correct
        RespEval=1;
        Resp=1;
    else
        RespEval=0;
        Resp=0;
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
%% Just-noticable-difference Function
% One up, one down for now. Change if needed
% Make sure to jitter the scr.background noise so that the longer interval is not obvious by noise alone.
% Maybe remove backfground noise?

function [JND_results,resp_eval,JND_UD]=JNDfunction(scr,time,stim,JNDprac)

% Parameters
comparison_length=2; % start value for comparison interval length in s
standard_length=0.7; % length of standrad interval (ideally the same as the interval shown in the main task)
init_adjustment_steps=0.3; % bigger steps in the beginning (first descent)
adjustment_steps=0.02; % size of adjustment each step in seconds
if ~stim.reducetrials % normal version
maxreversals=18;
else % reduced trial N when testing the code
    maxreversals=5;
end
avg_trials=5; % average across how many of the last trials to calculate JND?

% Initialize
jndfound=0;
JND_results=table(); %preallocate results
currtrial=0;

% Initialize Palamedes Staircase Procedure
up = 1;                     % inceasypracticerease after 1 wrong
down = 1;                   % decrease after 3 consecutive right
stopcriterion = 'reversals';
stoprule = maxreversals; % needs to be defined, but I don't want it to stop by itself so I set it to higher than the trial number
startvalue = comparison_length;           %intensity on first trial
xMin=standard_length; % minimum length

if ~JNDprac % if not practice, run until max reversals are done
    UD = PAL_AMUD_setupUD('up',up,'down',down,'StepSizeDown',init_adjustment_steps,'StepSizeUp', ...
        init_adjustment_steps,'stopcriterion',stopcriterion,'stoprule',stoprule, ...
        'startvalue',startvalue, 'xMin',xMin);
else % if practice, run only one trial
    UD = PAL_AMUD_setupUD('up',up,'down',down,'StepSizeDown',init_adjustment_steps,'StepSizeUp', ...
        init_adjustment_steps,'stopcriterion','trials','stoprule',1, ...
        'startvalue',startvalue, 'xMin',xMin);
end

% Shuffle textures
cuetex=stim.cuetex(randperm(length(stim.cuetex)));
noisetex=stim.noisetex(randperm(length(stim.noisetex)));

% % JND Start Trigger
io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{1,'JND Start'});
pause(0.1)
io64(scr.triggerPort, scr.triggerPortAddress,0);

while ~UD.stop

    currtrial=currtrial+1; %Update trial counter
    cnoise=1;
    ccue1=1;
    ccue2=1;

    % Shuffle which interval gets presented first
    interval_lengths=[round(standard_length/scr.ifi), round(UD.xCurrent/scr.ifi)]; % combine and convert into frames
    interval_lengths=interval_lengths(randperm(length(interval_lengths)));

    % Shuffle III and post-interval mask lengths
    time.JNDpostintmask=time.JNDpostintmask(randperm(length(time.JNDpostintmask)));
    time.III_JND=time.III_JND(randperm(length(time.III_JND)));

    % Show Intervals
    for intervals=1:length(interval_lengths)
        % Fixation
        for fixidx=1:time.JNDpostintmask(intervals)
            Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour); % Adds current interval number to screen
            Screen('Flip', scr.win);
            cnoise=cnoise+1;
                        if fixidx==1
                            io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{intervals,'JND Mask Onset'});
                        elseif fixidx==2
                            io64(scr.triggerPort, scr.triggerPortAddress,0);
                        end
        end

        % Interval
        for cueidx=1:time.cuedur
            Screen('DrawTexture', scr.win, cuetex(ccue1), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
            Screen('Flip', scr.win);
            ccue1=ccue1+1;
                        if cueidx==1
                            io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{intervals,'JND Cue 1'});
                        elseif cueidx==2
                            io64(scr.triggerPort, scr.triggerPortAddress,0);
                        end
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
                        if cueidx==1
                            io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{intervals,'JND Cue 2'});
                        elseif cueidx==2
                            io64(scr.triggerPort, scr.triggerPortAddress,0);
                        end
        end

        % Fixation
        for fixidx=1:time.JNDpostintmask(intervals+1)
            Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(scr.win,int2str(intervals),'center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
            Screen('Flip', scr.win);
            cnoise=cnoise+1;
        end

        % Blank screen to mark end of interval
        for fixidx=1:time.III_JND(intervals+1)
            Screen('Flip', scr.win);
        end
    end

    % Query Response
        io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{1,'JND Question'});

    [response]=respfunction(scr.win,'Which interval was longer? \n\n 1 or 2?',["1","2"]);

        io64(scr.triggerPort, scr.triggerPortAddress,0);

    % Evaluate Response
    if  interval_lengths(1)==round(UD.xCurrent/scr.ifi) && UD.xCurrent>standard_length  && double(response)==1 % if (comparison was presented first & is longer) AND (the answer was 1) --> correct
        resp_eval=1;
                io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{1,'JND Response 1 Correct'});
    elseif interval_lengths(2)==round(UD.xCurrent/scr.ifi) && UD.xCurrent>standard_length && double(response)==2 % if (comparison was presented second & is longer) AND (the answer was 2) --> correct
        resp_eval=1;
                io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{1,'JND Response 2 Correct'});
    elseif interval_lengths(2)==round(UD.xCurrent/scr.ifi) && UD.xCurrent==standard_length % standard and comparison are equal
        resp_eval=0; % always incorrect when comparison and standard are equal
    else
        resp_eval=0; % incorrect if none of the above hold
    end

        io64(scr.triggerPort, scr.triggerPortAddress,0);

    % Update staircase
    UD = PAL_AMUD_updateUD(UD, resp_eval);

    % Update step size after third reversal (to narrow the steps, assuming the participant has now reached a value around their threshold)
    if ismember(3,UD.reversal)
        UD = PAL_AMUD_setupUD(UD,'StepSizeDown',adjustment_steps,'StepSizeUp', ...
            adjustment_steps);
    end

    % Register data
    JND_results(currtrial,:)={standard_length, UD.xCurrent, response, resp_eval};
    JND_UD=UD;

end
JND_results.Properties.VariableNames={'Standard','Comparison','Response','Correct'}; % add variable names to results table

% % JND Stop Trigger
io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervectorJND{1,'JND Stop'});
pause(0.1)
io64(scr.triggerPort, scr.triggerPortAddress,0);
end
%% Timing reminder Function
% Reminder Intervals at full visibility
function timing_reminder(scr, time, stim, trialinfoin)

% Randomly choose from jittered timings
ITI=time.ITIreminder(randi(length(time.ITIreminder)));
preq=time.preq(randi(length(time.preq)));
initmask=time.initmask(randi(length(time.initmask)));

% Select regular timing of the current condition as ISI (first trial is always regular)
if trialinfoin(3)==1
    time.ISI=round(time.trialtimes1(1)/scr.ifi); % convert current target interval to frames
elseif trialinfoin(3)==2
    time.ISI=round(time.trialtimes1(2)/scr.ifi); % convert current target interval to frames
else
    error("Error in selecting the timing reminder interval")
end

% Shuffle textures and initialize counters
cnoise=1;
ccue=1;
cuetex=stim.colCueTex(randperm(length(stim.colCueTex)));
noisetex=stim.noisetex(randperm(length(stim.noisetex)));

%     % Text
%     DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
%     Screen('Flip', scr.win);
%     pause (1)

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
%         io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'Reminder Interval Mask Onset'});
%     elseif fixidx==2
%         io64(scr.triggerPort, scr.triggerPortAddress, 0);
%     end
    cnoise=cnoise+1;
end

% Cue Onset
for cueidx=1:time.cuedur
    DrawFormattedText(scr.win,'Remember this interval','center', scr.yCenter-(0.25*scr.axisy(2)),scr.fontcolour);
    Screen('DrawTexture', scr.win, cuetex(ccue), [], [], 0);
    Screen('Flip', scr.win);
%     if cueidx==1
%         io64(scr.triggerPort, scr.triggerPortAddress,stim.triggervector{trialinfoin(3),'Reminder Interval Cue 1'}); % cue trigger
%     elseif cueidx==2
%         io64(scr.triggerPort, scr.triggerPortAddress, 0);
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
%         io64(scr.triggerPort, stim.triggervector{trialinfoin(3),'Reminder Interval Cue 2'});
%     elseif taridx==2
%         io64(scr.triggerPort, scr.triggerPortAddress, 0);
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
