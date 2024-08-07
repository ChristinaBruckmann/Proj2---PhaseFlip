%% Plot figure in PTB
% Christina Buckmann 07.08.2024

% Takes a matlab figure and converts it into a texture, before plotting it to the current PTB window.

% Requires PTB to be initialized and a window to be open
% Requires a figure to already be plotted.

% Potential improvements: 
% Allow for position input, currently plots in the center
% Needs also a way to adjust figure size (might be combined with position input?)

function PTB_plotfig(figure_handle, screen_handle, save_name, flipnow)

% Argument Handling
arguments
    figure_handle % Matlab figure handle - which figure to plot?
    screen_handle % PTB screen handle - which open window to plot to?
    save_name string ="Figure_to_plot" % Optional, if not given replace by default
    flipnow {mustBeMember(flipnow, [0, 1])} =1 % Flip right away inside function? Must be 0 or 1, default is 1
end

% Save figure as png and close
saveas(figure_handle, save_name,'png');
close(figure_handle);

% Load the image file
load_name=strcat(save_name,'.png');
img = imread(load_name);

% Make the texture
texture = Screen('MakeTexture', screen_handle, img);

% Draw the texture on the screen
Screen('DrawTexture', screen_handle, texture);

if flipnow
    % Flip to the screen
    Screen('Flip', screen_handle);
    KbStrokeWait;
end

end
