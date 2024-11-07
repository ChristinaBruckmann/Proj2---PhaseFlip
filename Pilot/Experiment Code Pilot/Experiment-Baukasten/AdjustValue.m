%% Adjust Contrast Level Function

%% PTB Set Up

% Setup PTB
PsychDefaultSetup(2);
Priority(1);
Screen('Preference', 'SkipSyncTests', 1);
screeninfo.number = max(Screen('Screens'));

% Screen Values
white = WhiteIndex(screeninfo.number);
screeninfo.black = BlackIndex(screeninfo.number);
screeninfo.grey = white / 2;
background=screeninfo.grey;
screeninfo.fontcolour=screeninfo.black;

[screeninfo.win, screeninfo.windowRect] = PsychImaging('OpenWindow', screeninfo.number, background, [], 32, 2,...
    [], [],  kPsychNeed32BPCFloat);

% Initial Intensity
gaborcontrast=0.4;

% Run function
gaborcontrast=adjustintensity(screeninfo,gaborcontrast);

sca;
disp('finished')

%% Adjusting Function
% Insert current intensity, adjust and return new intensity.
% Increase/Decrease with arrows, confirm with enter
function [newintensity]=adjustintensity(screeninfo, gaborcontrast)
adjusting=1;
currlevel=gaborcontrast;
    while adjusting
    currentleveltext=sprintf('Current luminance level: %.2f',currlevel);
    DrawFormattedText(screeninfo.win,currentleveltext, 'center', 'center', 1);
    Screen('Flip', screeninfo.win);
    
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