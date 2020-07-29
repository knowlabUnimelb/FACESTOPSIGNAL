function closeexp (screenparms)
% CLOSEEXP Close Experiment Window
%   CLOSEEXP (screenparms) closes current window pointed to by
%   handles in screenparms
%
%   see also prepexp
%   created by Stephan Lewandowsky 2003

ShowCursor; % show Mouse Cursor
% % Screen(screenparms.window,'Close'); % close Experiment Window
Screen('CloseAll')
ListenChar(1)