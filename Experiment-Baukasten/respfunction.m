%% Response Query and Graceful Abort
% Christina Bruckmann 07.08.2024
% Still alpha-testing
% Planned improvements: optional RT output
%
% For Psychtoolbox
% PTB Window needs to be opened
%
% Input
    % Format: [response,resp_eval]=respfunction(screen_handle,'Test',["RightArrow","LeftArrow"],[0,1],["9","8"],0);
    % Required:
        % Screen Handle: From PTB
        % Question Text: as char. Accepts sprintf format (use single quotes not double, otherwise crashes)
        % Response Buttons: String array, as many options as needed, separated by commas. Must be at least two.
    % Optional
        % Response Mapping: Double-Vector, which responses are counted as correct (also other mappings possible. Output: resp_eval
        % Graceful Abort Buttons: Default e - escape and c-continue. 
        % Fontcolour im PTB format. 
%
% Output
        %Response: The button that was pressed (as a string, even when numbers are used as buttons)
        %resp_eval: Corresponding response mapping, should it have been provided. Otherwise returns empty.

function [response,resp_eval]=respfunction(screen_handle,questiontext,responsebuttons,responsemapping,graceful_abort,fontcolour)

%% Check input
arguments
    screen_handle
    questiontext char 
    responsebuttons (1,:) string =["1","2"]; % Response buttons
    responsemapping (1,:) double =0; % Response mapping - which buttons are correct(1)/incorrect(0) 
    graceful_abort (1,2) string =["e","c"]; % Escape and continue. default e and c
    fontcolour (1,1) double  = 0 % default black
end

% At least two response buttons declared?
 if ~(numel(responsebuttons)>1)
    error("Please provide at least two response options.")
end

% Buttons are valid keys?
allkeys=KbName('KeyNames'); % get all possible keys
allkeys=allkeys(~cellfun(@isempty, allkeys)); % remove empty fields

for nbuttons=1:length(responsebuttons)
    if ~(ismember(responsebuttons(nbuttons),allkeys))
        error ("Response option in position %i is invalid.\nFor a list of all options run KbName('KeyNames')",nbuttons)
    end
end

% No duplicates in response buttons
catArray = categorical(responsebuttons);% Convert the cell array of strings to a categorical array
[uniqueValues, ~, uniqueIdx] = unique(catArray); % Find the unique values and their counts
counts = accumarray(uniqueIdx, 1);
if any(counts > 1) % Check if any counts are greater than 1
    duplicateValues = uniqueValues(counts > 1);  % Find the duplicate values
    duplicateStr = strjoin(cellstr(duplicateValues), ', '); % Convert duplicate values to a comma-separated string for the error message
    error('Duplicate response options: %s. Each button has to be uniquely assigned.', duplicateStr);
end

% If response mapping is provided: Number of response buttons must match number of response mappings
if sum(responsemapping)>0 && ~(numel(responsebuttons)==numel(responsemapping))
    error ("Number of response buttons must match number of response mappings.")
end

% Graceful abort valid buttons
for nbuttons=1:length(graceful_abort)
    if ~(ismember(graceful_abort(nbuttons),allkeys))
        error ("Graceful abort input in position %i is invalid.\nFor a list of all options run KbName('KeyNames')",nbuttons)
    end
end

% Not overlap between buttons
if ~isempty(intersect(responsebuttons,graceful_abort))
    error ("Response buttons and graceful abort buttons cannot be the same key.")
end

%% Query response
resp=0; 

while~resp
    DrawFormattedText(screen_handle,questiontext, 'center', 'center',fontcolour);
    Screen('Flip', screen_handle);
    [~, keyNamethr, ~]=KbStrokeWait;

    % Graceful abort? 
    if strcmp(KbName(keyNamethr),graceful_abort(1))==1
        abort_text=sprintf('Continue(%s)? Exit(%s)?',graceful_abort(2),graceful_abort(1));
        DrawFormattedText(screen_handle,abort_text, 'center', 'center', fontcolour);
        Screen('Flip', screen_handle);
        [~, keyNamethr, ~]=KbStrokeWait;
        if strcmp(KbName(keyNamethr),graceful_abort(1))
            ListenChar(0);
            sca
            return
        elseif strcmp(KbName(keyNamethr),graceful_abort(2))
            continue
        end
    end

    % If not, check all possible response options
    for nbuttons=1:length(responsebuttons)
        if strcmp(KbName(keyNamethr),responsebuttons(nbuttons))==1
            response=responsebuttons(nbuttons); % register response
            % If response mappings are provided, also evaluate correctnes of response.
            if ~sum(responsemapping)==0
                resp_eval=responsemapping(nbuttons);
            else
                resp_eval=[];
            end
            resp=1;
        end
    end
end
end