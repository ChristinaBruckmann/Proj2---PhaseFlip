%% Test New Baukasten Functions

%% Initialize PTB
floatwin=1; % Take over whole screen (0) or just part of it (1)
SkipSync=1; % Skip Synch (1 when testing on laptop)

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

%% Test Function here
try
    DrawFormattedText(scr.win,'Do not change this slide.', 'center', 'center', scr.fontcolour);
    Screen('Flip', scr.win);
    unlock_continue(scr.win,["123"])
catch
    sca;
    ShowCursor;
    psychrethrow(psychlasterror);
end