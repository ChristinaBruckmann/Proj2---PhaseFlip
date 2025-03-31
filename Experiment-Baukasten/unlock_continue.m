%% Preventing Participants from skipping ahead by themselves
% Requires Specific Numer/Letter-Code Input to Continue
% Enter code and confirm with enter. If code is correct, function terminates and returns to script.
%
% Christina Bruckmann 13.09.2024
%
% For Psychtoolbox
% PTB Window needs to be opened

%%%%% IMPORTANT: GetString() tends to be buggy when mex files are cleared.%%%%%%%

%
% Input:
    % Required
            % screen_handle: screen handle from PTB
            % disp_text: string that is displayed on screen
            % unlock_code: string array that is required in that specific order for continuing the script
    % Optional
            % Graceful Abort Buttons: Default e - escape and c-continue. 
            % Fontcolour in PTB format, default black.
            % background colour in PTB format, default grey.


%% Function
function unlock_continue(screen_handle, unlock_code,disp_text,graceful_abort,fontcolour)

% Input Checks
arguments
    screen_handle
    unlock_code (1,:) string % Response buttons
    disp_text char =[]; % If display text is provided, the code prints text to screen. If no text is provided, it will not change the previous screen until it is unlocked
    graceful_abort (1,2) string =["e","c"]; % Escape and continue. default e and c
    fontcolour (1,1) double  = 0 % default black
end

% Graceful abort valid buttons
allkeys=KbName('KeyNames'); % get all possible keys
allkeys=allkeys(~cellfun(@isempty, allkeys)); % remove empty fields

for nbuttons=1:length(graceful_abort)
    if ~(ismember(graceful_abort(nbuttons),allkeys))
        error ("Graceful abort input in position %i is invalid.\nFor a list of all options run KbName('KeyNames')",nbuttons)
    end
end

% Not overlap between buttons
if ~isempty(intersect(unlock_code,graceful_abort))
    error ("Unlock code cannot include graceful abort buttons. \n\n Default graceful abort buttons are e and c.")
end

%% Run Function
unlocked=0;

while ~unlocked % As long as the correct code has not been entered:
    % Wait for correct code
    if ~isempty(disp_text) % if display text has been provided
    DrawFormattedText(screen_handle,disp_text, 'center', 'center', fontcolour);
    Screen('Flip', screen_handle);
    end
    entered_code = GetString();
    % Graceful abort?
    if strcmp(entered_code,graceful_abort(1))==1
        if ~isempty(disp_text)
        abort_text=sprintf('Continue(%s)? Exit(%s)?',graceful_abort(2),graceful_abort(1));
        DrawFormattedText(screen_handle,abort_text, 'center', 'center', fontcolour);
        Screen('Flip', screen_handle);
        end
        [~, keyNamethr, ~]=KbStrokeWait;
        if strcmp(keyNamethr,graceful_abort(1))
            sca
            return
        elseif strcmp(keyNamethr,graceful_abort(2))
            continue
        end
    end

    % If not, check if entered code matched unlock code. If not display a warning and restart loop.
    if strcmp(entered_code,unlock_code)==1
        unlocked=1; % unlock and continue with script
        Screen('Flip', screen_handle);
    else
        if ~isempty(disp_text)
            DrawFormattedText(screen_handle,'Incorrect code. \n\n Please ask the experimenter to continue. \n\n If you are the experimenter, press any button to re-enter the code.', 'center', 'center', fontcolour);
            Screen('Flip', screen_handle);
            KbStrokeWait;
        end
    end

end
end