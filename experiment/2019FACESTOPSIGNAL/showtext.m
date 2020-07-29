% SHOWTEXT Display text
%       To be used with Psychophysics Toolbox .mex file, SCREEN
%       SHOWTEXT (screenparms, size, stimulus, onoff, vertdown, leftright)
%       displays 'stimulus' in a window defined in 'screenparms' where
%       'size' = stimulus size in pixels
%
%       'onoff' specifies color (0 = white, 1 = black)
%       'vertdown' specifies vertical offset from center in pixels
%       'leftright' specifies horizontal offset from center in pixels
%
%       See also prepexp, screen



function showtext (screenparms, size, stimulus, onoff, vertdown, leftright)
global ptb3

size = round(size);
if onoff
    c = screenparms.black;
else
    c = screenparms.white;
end

% fromtop = (screenparms.rect(RectBottom) - screenparms.rect(RectTop) + size)/2 + vertdown;
fromtop = (screenparms.rect(RectBottom) - screenparms.rect(RectTop))/2 + vertdown;
fromleft = (screenparms.rect(RectRight) - screenparms.rect(RectLeft))/2 - (length(stimulus)/2 * size/2.4) + leftright;

if ptb3
    Screen('TextSize', screenparms.window,size);
    Screen('TextFont', screenparms.window,screenparms.sansSerifFont);
    Screen('DrawText', screenparms.window, stimulus, fromleft, fromtop, c);
else
    Screen(screenparms.window,'TextSize',size);
    Screen(screenparms.window,'TextFont',screenparms.sansSerifFont);
    Screen(screenparms.window,'DrawText', stimulus, fromleft, fromtop, c);
end