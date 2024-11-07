%% Adjust Value Function
% Insert current value, adjust and return new value.
% Increase/Decrease with arrows, confirm with enter
function [newvalue]=adjustintensity(screen_handle, oldvalue,fontcolour)
adjusting=1;
currlevel=oldvalue;
    while adjusting
    currentleveltext=sprintf('Current value: %.2f',currlevel);
    DrawFormattedText(screen_handle,currentleveltext, 'center', 'center', fontcolour);
    Screen('Flip', screen_handle);
    
    [keyDown, ~, keyName] = KbCheck;
    if keyDown
        if strcmp(KbName(keyName), 'UpArrow')==1
            currlevel=currlevel+0.01;
            pause(0.05)
        elseif strcmp(KbName(keyName), 'DownArrow')==1
            currlevel=currlevel-0.01;
            pause(0.05)
        elseif strcmp(KbName(keyName), 'Return')==1 % confirm with enter
            newvalue=currlevel;
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