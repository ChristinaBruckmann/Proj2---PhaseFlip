%% Presents a GUI menu for the user to interact with
% First input is the title of the menu, then a variable amount of choice options can be passed
% Returns the value given by the user, if the window is closed it returns 0
% Written by ChatGPT :)

function choice = centeredMenu(title, varargin)
    % Check if at least one option is provided
    if nargin < 2
        error('You must provide at least one option.');
    end

    % Screen size
    screenSize = get(0, 'ScreenSize');
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);

    % Determine the longest string for resizing the dialog width
    allStrings = [{title}, varargin];  % Include title in the string length calculation
    maxStringLength = max(cellfun(@length, allStrings));

    % Set basic dimensions
    minDialogWidth = 300;  % Minimum dialog width
    charWidthFactor = 8;   % Estimated character width in pixels
    dialogWidth = max(minDialogWidth, charWidthFactor * maxStringLength);  % Adjust width

    % Dialog height depends on the number of options
    numOptions = length(varargin);
    baseHeight = 100;
    optionHeight = 50;
    dialogHeight = baseHeight + numOptions * optionHeight;

    % Calculate position for centering the window
    dialogPosX = (screenWidth - dialogWidth) / 2;
    dialogPosY = (screenHeight - dialogHeight) / 2;

    % Create figure for the dialog box with dynamic width and height
    hFig = figure('Position', [dialogPosX, dialogPosY, dialogWidth, dialogHeight], ...
        'MenuBar', 'none', 'Name', title, 'NumberTitle', 'off', ...
        'Resize', 'off', 'WindowStyle', 'modal', 'CloseRequestFcn', @cancelCallback);

    % Initialize choice
    choice = 0;

    % Create title text
    uicontrol('Style', 'text', 'String', title, ...
        'Position', [50, dialogHeight - 80, dialogWidth - 80, 40], ...
        'HorizontalAlignment', 'center', 'FontSize', 12);

    % Create buttons dynamically based on the number of options
    for i = 1:numOptions
        uicontrol('Style', 'pushbutton', 'String', varargin{i}, ...
            'Position', [50, dialogHeight - 90 - i*optionHeight, dialogWidth - 100, 40], ...
            'Callback', @(src,~) buttonCallback(i));
    end

    % Wait for the figure to close before returning the choice
    uiwait(hFig);

    % Callback function for options
    function buttonCallback(optionIdx)
        choice = optionIdx;
        delete(hFig); % Close dialog
    end

    % Callback for canceling the dialog (i.e., if the window is closed manually)
    function cancelCallback(~, ~)
        choice = 0; % Return 0 if the user closes the window without selecting an option
        delete(hFig); % Close the dialog
    end
end