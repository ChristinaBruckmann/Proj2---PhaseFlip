%% Instructions and Page Navigation
% Christina Bruckmann 06.08.2024

% Needs to be embedded in a code which has already PTB initialized and window opened etc.
% Allows to create instructions with back and forth navigations through button presses
% If desired also add position input to code

% Input: 
% screen: PTB screen handle
% instructions: string arrayformatted instructions, each new page should be a new row in the cell. amount of pages will be determined by row of cell array.
% optional: String array: buttons (which buttons should mean forward and backward? First is forward, second is backward) Default: arrow left right e.g. ["f","b"]
% not yet functional: Optional: page counter at the bottom (e.g. 4/6) - yes(1), no(0) default no


function navipage(screen_win,instructions,fontcolour,buttons)

%% Argument handling

arguments
    % Required
    screen_win
    instructions (:,:) cell % need to be a cell array or character vectors (use single quotes not double, otherwise crashes)
    % Optional
    fontcolour (1,1) double  = 0
    buttons (1,:) string =["LeftArrow","RightArrow"];  % needs to be a string array with one row and as many characters as needed. Default buttons are the left and right arrow keys.
    %page_counter {mustBeMember(page_counter, [0, 1])}  =0% Must be 0 or 1 (default 0)
end

% Two buttons declared?
 if ~(numel(buttons)==2)
    error("Number of buttons must be exactly two.")
end

% Buttons are valid keys?
allkeys=KbName('KeyNames'); % get all possible keys
allkeys=allkeys(~cellfun(@isempty, allkeys)); % remove empty fields

if ~(ismember(buttons(1),allkeys) && ismember(buttons(2),allkeys))
    error ("Buttons must be valid keys.\nThey must be entered as a string array and match the keynames from Psychtoolbox.\nFor a list of all options run KbName('KeyNames')")
end

%% Run function

curr_slide=1;
final_slide=0;
while ~final_slide
    resp=0; % responded?
    while~resp % Wait for button press

        % Show current instruction slide
        DrawFormattedText(screen_win,instructions{curr_slide}, 'center', 'center', fontcolour);
        Screen('Flip', screen_win);

        % Wait for response
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
        elseif strcmp(KbName(keyNamethr),buttons(1))==1
            if ~curr_slide==1
            curr_slide=curr_slide-1; % previous page
            end
            resp=1;
        elseif strcmp(KbName(keyNamethr),buttons(2))==1 
            curr_slide=curr_slide+1; % next page
            resp=1;
        end
    end

    % Exit after last page
    if curr_slide>size(instructions)
        final_slide=1;
    end
end
end