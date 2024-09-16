%% Experiment Code Project 2 - Phase Flip
% To do:

        
    % Minor Fixes:
        % jitter duration for noise presentation before and after each trial
        % Add file name 
        % subresults saves timing as frames, change to seconds

    % Major Steps:
        % Add Instructions
        % Add practice!! (for practice also add "practice label "
        % Include GUI and Option to re-start in the middle.
        % Test code functionality for the whole experiment
        % EEG Prep Code

    % Discuss with Assaf:
        % Constraints for Matrix Shuffle
        % Review timing of different trial types (random values used so far)
         % Add gamification?


% Requires the custom functions  PTB_plotfig.m and navipage.m, unlock_continue.m and respfunction.m

% Clear
sca;
close all;
clear;
%clear mex;

% Initialize Structs
scr=[]; % Everything related to PTB Screen
stim=[]; % Stimulus Information
time=[]; % Timing Information

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj2 - PhaseFlip\Pilot\Experiment Code Pilot'
%% Input GUI
input_complete=0;
while ~ input_complete
    % Which type of settings would you like?
    title = 'What would you like to do?';
    options = {'Run Participant', 'Trigger Check', 'Oscilloscope Test','Debugging'};

    run_type = centeredMenu(title, options{:});

    switch run_type
        case 0
            f=errordlg("Please select what you would like to do.");
            uiwait(f)
        case 1 % Running a participant
            floatwin=0; % Take over whole screen (0) or just part of it (1)
            SkipSync=1; % Skip Sync Tests (1 when testing on laptop)
            testing=0; % testing the code? reduces amount of trials to a minimum
            speedrun=0; % automatically chooses a response, no need to manually click anything (used for testing the code)
            scr.scope=0; % run scope test?

            % Experiment Sections
            JND_task=1; % JND Task?
            staircase=1; % Staircase?
            practice=1; % Practice?

            % Input Subject Number
            sub_num_compl=0;
            subresults.subj_num=[];
            valid_subnumber=0;
            while ~sub_num_compl
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
                    end
                end

                % Check if a file for this participant already exists
                savefilename=sprintf("Pilot_PhaseFlip_Subj%i",subresults.subj_num);
                try
                    load(savefilename) % Try to load a file with this name, if it exists ask if this should be used
                    answer = questdlg('Subject already exists. Load existing file?','Existing Subject', 'Yes','Return.','Return.');
                    switch answer
                        case 'Yes'
                            % Start where the participant has left off?
                            answer = questdlg('Reloading file. Start where the participant left off?','Existing Subject', 'Yes','No - restart.');
                            switch answer
                                case 'Yes'
                                     sub_num_compl=1;
                                case 'No - restart.'
                                    old_subresults=subresults; % save previous results in file under new name
                                    save(savefilename,old_subresults)
                                    subresults=[]; % restart with clean subresults
                                    save(savefilename,subresults)
                                    f=msgbox("New file created.");
                                    uiwait(f)
                                    sub_num_compl=1;
                            end
                        case 'Return'
                            continue % Re-enter subject number
                    end
                catch ME
                    answer = questdlg('New subject. Create a new file?','New Subject', 'Yes','Return','Return');
                    % Handle response
                    switch answer
                        case 'Yes'
                            save(savefilename,"subresults");
                            f=msgbox("New file created.");
                            uiwait(f)
                            sub_num_compl=1;
                        case 'Return'
                            continue % Re-enter subject number
                    end
                end
            end

            % Input Experimenter Name
            subresults.experimenter = {};
            while true
                answer = inputdlg('Experimenter Name:');
                if isempty(answer)  % User closed the dialog without entering anything
                    continue;       % Reopen the dialog
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
            testing=0;
            speedrun=0;
            scr.scope=0;

            % Experiment Sections
            JND_task=0; 
            staircase=0; 
            practice=0; 

            savefilename=sprintf('trigger_test_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            input_complete=1;
        case 3 % Oscilloscope Test
            floatwin=0; 
            SkipSync=0; 
            testing=0;
            speedrun=0; 
            scr.scope=1;

            input_complete=1;
        case 4 % Debugging
            % Choose which debugging options you want
            prompt = 'Select debugging options:';
            options = {'Reduced_TrialN', 'Speedrun', 'JNDTask','Staircase','Practice','SkipSync','Floating_Window'};

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
            
            savefilename=sprintf('test_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            scr.scope=0;
            input_complete=1;
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
scr.background=scr.grey;
scr.fontcolour=scr.black;

% Open Screen (with or without floating window)
if ~floatwin
[scr.win, scr.windowRect] = PsychImaging('OpenWindow', scr.number, scr.background, [], 32, 2,...
    [], [],  kPsychNeed32BPCFloat);
else
    [scr.win, scr.windowRect] = PsychImaging('OpenWindow', scr.number, scr.background,...
    [100 100 1600 1600], [], [], [], [], [], kPsychGUIWindow);
end
scr.ifi = Screen('GetFlipInterval', scr.win);

% Define Screen Dimensions
[scr.xCenter, scr.yCenter] = RectCenter(scr.windowRect);
[scr.axisx]=scr.windowRect([1,3]);
[scr.axisy]=scr.windowRect([2,4]);

%% Timing Parameters
% Target interval parameters can be found in the trial matrix creation

time.ITI=1; % Inter-Trial Interval
time.preq=0.5; % Pre-Question Interval
time.cuedur=0.1; % Cue Duration
time.tardur=0.016; % Target Duration
     
% Target times for each block type / condition
time.trialtimes1=[0.7, 0.75, 1.2, 0.5, NaN]; %[timereg timeirr timelong timeshort timecatch];
time.trialtimes2=[0.75, 0.7, 1.25, 0.55, NaN]; %[timereg timeirr timelong timeshort timecatch];

% Convert into frames
time.ITI=round(time.ITI/scr.ifi); % Inter-Trial Interval
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
gaborpercent=1-stim.maskintensity; % Initialize like this, adjust with staircase later

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

text.JND_instructions={'JND_Task. \n\n Continue with arrows.';'Instructions Page 1 \n\n Continue with arrows.';'Instructions Page 2 \n\n Continue with arrows.';'Instructions Page 3 \n\n Continue with arrows.'};


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
    cond=[1 2 1 2 1 2 1 2]; % conditions (1=800ms, 2=850ms)
    nblocks=length(cond); % total blocks
    ntrials=50; % trials per block
    ntrialstot=nblocks*ntrials;

    % Split per block
    nregular=35;
    nirregular=3;
    nlong=1;
    nshort=1;

    % Left overs are assigned to catch trials (throw error if none are left)
    if (ntrials-nregular-nirregular-nlong-nshort)>49
        error('No catch trials left. Double check trial matrix generation. Catch trials are the most important part of this project.')
    elseif (ntrials-nregular-nirregular-nlong-nshort)>45
        warning('Using less than 5 catch trials per block. This is not advised. Please reconsider.')
    else
        ncatchtrials=ntrials-nregular-nirregular-nlong-nshort;
    end
end

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
%% Create Stimuli

maxframesnoise=time.ITI+time.ITI+time.preq; % How many frames are needed per trial
maxframescue=time.cuedur;
DrawFormattedText(scr.win, 'Loading. Please wait.', 'center', 'center', scr.fontcolour);
Screen('Flip', scr.win);

%  Noise Mask
for fridx=1:maxframesnoise*4 % make 4 times the amount of frames just to be sure and also for JND task more are needed.
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
 % Embedding it in mask happens lower after the intensity has been chosen
 

 % %% Oscilloscope Test
 % DrawFormattedText(scr.win, 'Welcome to the Oscilloscope Testing Environment \n\n Press any button to continue.', 'center', 'center', scr.fontcolour);
 % Screen('Flip', scr.win);
 % responsetext=sprintf("Which condition would you like to test? \n\n\n 1) %i seconds \n\n 2) %i seconds ",)
 % [response,resp_eval]=respfunction(scr.win,responsetext,responsebuttons,responsemapping,graceful_abort,fontcolour)
 %% Save all parameters and other info
 subresults.trialmatrix=trialmatrix;
 subresults.screeninfo=scr;
 subresults.textinfo=text;
 save(savefilename, 'subresults')

 %% Oscilloscope Test
 if scr.scope
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
     trialfunction(scr,time,text,stim,gaborpercent,[99,scope_dur],0); % Run trial

     end
 end
%% Run Experiment
try

    % Suppress key presses in command window
    ListenChar(2);

     % Introduction
    DrawFormattedText(scr.win,'Welcome to the experiment. \n\n Press any button to start.', 'center', 'center', scr.fontcolour);
    Screen('Flip', scr.win);
    KbStrokeWait;

    % JND Task
    while JND_task
        [JND, JND_results]=JNDfunction(scr,time,text.JND_instructions,stim);
        % Plot results and flip to screen
        Screen('Flip', scr.win);
        JND_fig = figure('Visible', 'off'); % make invisible figure
        plot(1:length(JND_results.Comparison), JND_results.Comparison, '-o', 'LineWidth', 2); xlim([1 length(JND_results.Comparison)]); ylim([min(JND_results.Comparison)-1 max(JND_results.Comparison)+1]); title('JND Results'); xlabel('Trial Number'); ylabel('Comparison Duration'); yline(JND) % plot data
        PTB_plotfig(JND_fig, scr.win, "JND_Figure", 0) % plot figure to PTB screen with this custom function
        Screen('Flip', scr.win);
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
            save(savefilename, 'subresults')
            JND_task=0; % Exit JND Task
        end
    end

    % Staircase
    if staircase
        [threshres,gaborpercent]=staircasefun(3,scr,time,text,stim,gaborpercent);
    end

    % Accept Staircase Output or Adjust Gabor Difficulty
    adjusttext=sprintf('Post-staircase intensity: %f \n\n Confirm (2) or Adjust (3)?',gaborpercent);
    [response]=respfunction(scr.win,adjusttext,["2","3"]);
    if double(response)==3 % if adjustment requested, ask for confirmation
        confirmationtext='Are you sure? \n\nIt is reccommended to stick to the staircase value. \n\n\n\n Accept staircase value (2) or Adjust anyway (3)?';
        [response]=respfunction(scr.win,confirmationtext);
        if double(response)==3
             gaborpercent=adjustintensity(scr.win, gaborpercent);
        end  
    end

    % Display the setting that will be used
    intconfirmed=sprintf('Confirmed Intensity: %f \n\n\n\nPress any button to continue.',gaborpercent);
    DrawFormattedText(scr.win,intconfirmed, 'center', 'center', scr.fontcolour);
    Screen('Flip', scr.win);
    KbStrokeWait;

    % Run Blocks
    tottrialcount=1; % total trial counter
  
    for b=1:nblocks

        % Determine current condition
        currcond=cond(b);

        % Run Trials
        for t=1:ntrials
            trialinfo=trialmatrix(tottrialcount,:); % Choose trialinfo from trialmatrix

            [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,gaborpercent,trialinfo,speedrun); % Run trial

            if ~isnan(Resp)% Save trial results
                resulttable(tottrialcount,:)=table(currcond, b, t, trialinfo(1), trialinfo(2),RT, Resp, RespEval, warning, gaborpercent,stim.maskintensity, 'VariableNames',{'Condition','Block', ...
                    'Trial','Trial Type','Target Interval','Reaction Time', 'Orientation Reseponse', 'Correct/Incorrect', 'Late Warning', 'Gabor Strength','Mask Intensity'});
                subresults.data=resulttable;
                save('test.mat', 'subresults')
                tottrialcount=tottrialcount+1; % update total trial counter
            else % repeat the trial and do not save output
                DrawFormattedText(scr.win,'Press any button to continue the task.', 'center', 'center', scr.fontcolour);
                Screen('Flip', scr.win);
                KbStrokeWait;
            end
        end

        % End of block/experiment message
        if b==nblocks
            DrawFormattedText(scr.win,'This is the end of the experiment. \n\n Please wait for the experimenter.', 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            unlock_continue(scr.win, scr.unlock_code) % Blocks screen until experimenter unlocks (to prevent subject from changing the slide)

            DrawFormattedText(scr.win,'Data has been saved. \n\n Press any button to close the screen.', 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            KbStrokeWait;
        else
            next_block=0; % continue with next block?
            while ~next_block
                blockendmessage=sprintf('End of block %i/%i \n\nPlease take a break. \n\n \n\nPress the space bar to continue.',b,nblocks);
                DrawFormattedText(scr.win,blockendmessage, 'center', 'center', scr.fontcolour);
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
%% Trial Function
function [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,gaborpercentin,trialinfoin,speed)
% trialinfoin is the input for the trial information(type of trial and target interval duration from the trial matrix)
% Choose Inter stimulus interval based on trialinfo input and convert to frames
% speedrun automatically chooses a response so no manual response is needed (for code testing)
time.ISI=round(trialinfoin(2)/scr.ifi);

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

% Start Trial
for fixidx=1:time.ITI
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    %Screen('DrawLines', scr.win, stim.fixCoords,stim.lineWidthPix, scr.black, [scr.xCenter scr.yCenter], 2);
    Screen('Flip', scr.win);
    cnoise=cnoise+1;
end

for cueidx=1:time.cuedur
    Screen('DrawTexture', scr.win, cuetex(ccue), [], [], 0);
    Screen('Flip', scr.win);
    ccue=ccue+1;
end

for noiseidx=1:time.ISI
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
    cnoise=cnoise+1;
end

for taridx=1:time.tardur
    Screen('DrawTexture', scr.win,  tartex(taridx), [], [], 0);
    Screen('Flip', scr.win);
end

for cueidx=1:time.preq
    Screen('DrawTexture', scr.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', scr.win);
    cnoise=cnoise+1;
end

% Request Orientation Reponse
startTime=GetSecs;
responded=0;
currtime=0;
warning=0;
if ~speed && ~scr.scope % if not in speed run mode or scope test, ask question as normal, otherwise automatically choose answer and continue
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
                else
                    RespEval=0;
                end
            elseif strcmp(KbName(keyName), text.keyright)==1 % if right key
                RT=respTime-startTime;
                responded=1;
                Resp=0; %Right
                if targetorient==0 %Right
                    RespEval=1;
                else
                    RespEval=0;
                end
            elseif strcmp(KbName(keyName),'e')==1
                DrawFormattedText(scr.win,'Continue (c)? Exit (e)?', 'center', 'center', scr.fontcolour);
                Screen('Flip', scr.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'e')
                    sca
                    return
                elseif strcmp(KbName(keyNamethr),'c')
                    Resp=NaN; % if trial was paused in between, do not go back to record late answers but continue with next trial instead
                    RespEval=0;
                    RT=NaN;
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
function [thresholdres,newgaborpercent]=staircasefun(reversals,scr,time,text,stim,oldgaborpercent)
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