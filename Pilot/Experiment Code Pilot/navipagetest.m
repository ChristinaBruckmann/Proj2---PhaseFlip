function navipagetest(screen, instructions, buttons)

arguments
    % Required
    screen
    instructions (:,:) string % Need to be a cell array with one column and as many rows as needed
    % Optional
    buttons (1,:) string = ["left", "right"];  % Needs to be a string array with one row and as many characters as needed. Default buttons are the left and right arrow keys.
end

% Debugging print statements
disp("Value of buttons: "), disp(buttons)

% Additional function implementation can go here

end