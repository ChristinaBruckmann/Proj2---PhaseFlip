%% Experiment Code Project 2 - Phase Flip
% Clear
sca;
close all;
clear;
clear mex;

% Initialize Structs
scr=[]; % Everything related to PTB Screen
stim=[]; % Stimulus Information
time=[]; % Timing Information
%% Usability Parameters
floatwin=1; % Take over whole screen (0) or just part of it (1)
scr.scope=0; % Oscilloscope Test?
SkipSync=1; % Skip Synch (1 when testing on laptop)
staircase=1; % Staircase?

% Blocks Trials etc.
cond=[1 2 1 2]; % conditions (1=800ms, 2=850ms)
nblocks=length(cond); % total blocks
ntrials=3; % trials per block
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
background=scr.grey;
scr.fontcolour=scr.black;

% Open Screen (with or without floating window)
if ~floatwin
[scr.win, scr.windowRect] = PsychImaging('OpenWindow', scr.number, background, [], 32, 2,...
    [], [],  kPsychNeed32BPCFloat);
else
    [scr.win, scr.windowRect] = PsychImaging('OpenWindow', scr.number, background,...
    [100 100 1600 1600], [], [], [], [], [], kPsychGUIWindow);
end
scr.ifi = Screen('GetFlipInterval', scr.win);
% Define Screen Dimensions
[scr.xCenter, scr.yCenter] = RectCenter(scr.windowRect);
[scr.axisx]=scr.windowRect([1,3]);
[scr.axisy]=scr.windowRect([2,4]);

%% Timing Parameters (in s)
time.ISI=0.8; % Target Interval
time.ITI=1; % Inter-Trial Interval
time.preq=0.5; % Pre-Question Interval
time.cuedur=0.1; % Cue Duration
time.tardur=0.016; % Target Duration

% Convert into frames
time.ISI=round(time.ISI/scr.ifi); % Target Interval
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

%% Create Stimuli

maxframesnoise=time.ISI+time.ITI+time.preq; % How many frames are needed per trial
maxframescue=time.cuedur;
DrawFormattedText(scr.win, 'Loading. Please wait.', 'center', 'center', scr.fontcolour);
Screen('Flip', scr.win);

%  Noise Mask
for fridx=1:maxframesnoise
    if scr.scope % If oscillosope test, make mask as dark as possible
        stim.noiseimg=0.01+stim.maskintensity*(rand(stim.rectSize)-0.5);
    else
        stim.noiseimg=0.5+stim.maskintensity*(rand(stim.rectSize)-0.5);
    end

    stim.noisetex(fridx)=Screen('MakeTexture', scr.win, stim.noiseimg); % Convert into Textures
end

% Cue Stimulus
[x, y]=meshgrid(-1:2/(stim.rectSize-1):1,-1:2/(stim.rectSize-1):1);
for fridx=1:maxframescue
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

%% Run Experiment
try

     % Introduction

    DrawFormattedText(scr.win,'Welcome to the experiment. \n\n Press any button to start.', 'center', 'center', scr.fontcolour);
    Screen('Flip', scr.win);
    KbStrokeWait;
    % Staircase
   [threshres,gaborpercent]=staircasefun(3,scr,time,text,stim,gaborpercent);

    % Accept Staircase Output or Adjust Gabor Difficulty
    adjusttext=sprintf('Post-staircase intensity: %f \n\n Confirm (2) or Adjust (3)?',gaborpercent);
    [response]=respfunction(adjusttext,scr);
    if response==0 % if adjustment requested, ask for confirmation
        confirmationtext='Are you sure? \n\nIt is reccommended to stick to the staircase value. \n\n\n\n Accept staircase value (2) or Adjust anyway (3)?';
        [response]=respfunction(confirmationtext,scr);
        if response==0
             gaborpercent=adjustintensity(scr, gaborpercent);
        end
    end

    % Display the setting that will be used
    intconfirmed=sprintf('Confirmed Intensity: %f \n\n\n\nPress any button to continue.',gaborpercent);
    DrawFormattedText(scr.win,intconfirmed, 'center', 'center', scr.fontcolour);
    Screen('Flip', scr.win);
    KbStrokeWait;

    % Run Blocks
    tottiralcount=1; % total trial counter
  
    for b=1:nblocks

        % Determine current condition
        currcond=cond(b);

        % Run Trials
        for t=1:ntrials
                [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,gaborpercent);

            % Save trial results
            resulttable(tottiralcount,:)=table(currcond, b, t, RT, Resp, RespEval, warning, gaborpercent,stim.maskintensity, 'VariableNames',{'Condition','Block', ...
                'Trial','Reaction Time', 'Orientation Reseponse', 'Correct/Incorrect', 'Late Warning', 'Gabor Strength','Mask Intensity'});
             subresults.data=resulttable;
             save(['test.mat'], 'subresults')
             tottiralcount=tottiralcount+1; % update total trial counter
        end

        % End of block/experiment message
        if b==nblocks
            DrawFormattedText(scr.win,'This is the end of the experiment. \n\n Please wait for the experimenter.', 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            KbStrokeWait;
        else
            blockendmessage=sprintf('End of block %i/%i \n\nPlease take a break. \n\n \n\nPress any button to start the next block.',b,nblocks);
            DrawFormattedText(scr.win,blockendmessage, 'center', 'center', scr.fontcolour);
            Screen('Flip', scr.win);
            KbStrokeWait;
        end
    end

    sca;
    ShowCursor;
catch
    sca;
    ShowCursor;
    psychrethrow(psychlasterror);
end
%% Trial Function
function [Resp, RespEval, RT, warning]=trialfunction(scr,time,text,stim,gaborpercentin)

    % Create gabors with the chosen Difficulty
    for i=1:time.tardur
        if scr.scope
            noiseimg=0.01+maskintensity*(rand(stim.rectSize)-0.5); % Dark noise mask
            % make cue area of noise tex white instead of showing Gabor
            xfrom=(length(stim.noiseimg)-stim.baseRect(3))/2;
            xto=xfrom+stim.baseRect(3);
            noiseimg(xfrom:xto,xfrom:xto)=stim.rectColorCue;
            stim.tartex(i)=Screen('MakeTexture', scr.win, noiseimg);
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
        tartex=stim.tartex;
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
            [~, respeval, ~, ~]=trialfunction(scr,time,text,stim,currentgabor);
            trialcount=trialcount+1;

            % Save
            thresholdres(threshrun,trialcount)=currentgabor;

            % Update
            if respeval==1 % if correct, increase streak
                streak=streak+1;
                if streak==3 % if third one in a row correct
                    currentgabor=currentgabor-stepsize; % decrease by step size
                    if lastchange=='inc'
                    reverscount=reverscount+1; % update reversal counter
                    end
                    lastchange='dec'; % register last change as decrease
                    if firstrev==0 % If this was the first reversal, count the trial on which this happened
                        firstrev=trialcount;
                    end
                    streak=0; % restart streak counter
                    prevstreak=1; % register that this was a streak
                end
            else % otherwise start from 0 and make task easier
                streak=0;
                prevstreak=0; % reset previous streak variable
                currentgabor=currentgabor+stepsize;
                if lastchange=='dec'
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
        end
    
        % Repeat?
        [response]=respfunction('Repeat threshold? (y-2/n-3)',scr);
        if response==1 % if repeat
            threshrun=threshrun+1;
        else % end staircase and save result
            finishstair=1;
        end

% Plot
        figure('Visible', 'off');
        plot(thresholdres(threshrun,:))
        yline(newgaborpercent)
        ylim([0 1])
        stairim = getframe; % Save figure as image
        stairtex=Screen('MakeTexture', scr.win, squeeze(stairim.cdata(:,:,1))); % Turn into texture
        Screen('DrawTexture', scr.win, stairtex, [], [], 0);
        Screen('Flip', scr.win);
        KbStrokeWait;
    end
end

%% Just-noticable-difference Function
% One up, one down for now. Change if needed

function [JND]=JNDfunction(scr,time,text,stim)

% Parameters
comparison_length=2; % start value for comparison interval length in s
standard_length=0.8; % length of standrad interval (ideally the same as the interval shown in the main task)
adjustment_steps=0.05; % size of adjustment each step in seconds
maxreversals=20;

% Initialize
jndfound=0;
currstep='none'; % preallocate previous response for reversal counter
revers_count=0; %preallocate reversal counter

while ~jndfound

    % Shuffle which interval gets presented first
    interval_lengths=[standard_length, comparison_length];
    interval_lengths=interval_lengths(randperm(length(interval_lengths)));

    % Show First Interval

    % Show Second Interval

    % Query Response

    % Evaluate Response
    if  interval_lengths(1)==comparison_length && comparison_length>standard_length  && response==1 % comparison is longer and was presented first and the answer was 1 --> correct
        resp_eval=1;
    elseif interval_lengths(2)==comparison_length && comparison_length>standard_length && response==2 % comparison is longer and was presented second and the answer was 2 --> correct
        resp_eval=1;
    elseif interval_lengths(2)==comparison_length && comparison_length==standard_length
        resp_eval=0; % always incorrect when comparison and standard are equal
    else
        resp_eval=0; % incorrect if none of the above hold
    end

    % Register data
    

    % Adjust Interval Length
    if  resp_eval==1 % response is correct
        % Lower comparison
        comparison_length_new=comparison_length-adjustment_steps; % lower comparison toward standard
        currstep='down'; % register which direction the current step was taken
        if comparison_length_new<standard_length % if the new comparison would be smaller than the standard (can happen depending on which step size is chosen, so that it 'jumps over' the equal stage), make equal to standard
            comparison_length_new=standard_length;
        end
    else % if response is incorrect, make comparison longer
        comparison_length_new=comparison_length+adjustment_steps;
        currstep='up'; % register which direction the current step was taken
    end

    % Reversal?
    if currstep~=prevstep % if current direction is not equal to previous direction
        revers_count=revers_count+1; % count as reversal
    end
    prevstep=currstep; % update step direction for next iteration

    % Update Comparison or JND Found? 
    if revers_count==maxreversals % if max reversals has been reached, take the average across the last 10 reversals as the JND
        jndfound=1;
        JND=;
    elseif corr_at_lowest_level==lim_lowest_level % if the answer at the lowest level was correct x times in a row, finish and register this as JND
        jndfound=1;
        JND=[standard_length comparison_length];
    else % if not found yet, update comparison
        comparison_length=comparison_length_new;
    end
end

end
%% Response Function
function [response]=respfunction(questiontext,scr)

resp=0; % responded?
        while~resp
            DrawFormattedText(scr.win,questiontext, 'center', 'center',scr.fontcolour);
            Screen('Flip', scr.win);
            [~, keyNamethr, ~]=KbStrokeWait;
            if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                DrawFormattedText(scr.win,'Continue (8)? Exit (9)?', 'center', 'center', scr.fontcolour);
                Screen('Flip', scr.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                    sca
                    return
                elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                    continue
                end
            elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                response=0;
                resp=1; 
            elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                response=1;
                resp=1; 
            end
        end

end