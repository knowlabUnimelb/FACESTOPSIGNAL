% Set up experiment using Psychophysics toolbox by opening on an onscreen
% window and setting fonts
function screenparms = prepexp(varargin)

warning off MATLAB:DeprecatedLogicalAPI
warning off MATLAB:mex:deprecatedExtension

%% Optional arguments
optargs = {0,...               % Screen Number {0 - full Windows desktop area, 1 - DISPLAY1, 2 - DISPLAY2
           [255, 255, 255],... % Color
           [],...              % Rect
           };
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[screenNumber, color, rect] = optargs{:}; % Place optional args in memorable variable names

%% Open window
screenparms = struct('window', [], 'rect', rect, 'color', color,...
    'x0', [], 'y0',[], 'topLeftQuad', [], 'topRightQuad', [],...
    'bottomLeftQuad', [], 'bottomRightQuad', [],...
    'white', [], 'black', [], 'serifFont', [], 'sansSerifFont', [],...
    'symbolFont', [], 'displayFont', [], 'screen', color);
[screenparms.window, screenparms.rect] =...
    Screen('OpenWindow', screenNumber, screenparms.color, rect);
[screenparms.x0, screenparms.y0] = RectCenter(screenparms.rect);

% Create screen quadrants
screenparms.topLeftQuad     = ScaleRect(screenparms.rect, .5, .5);
screenparms.topRightQuad    = OffsetRect(screenparms.topLeftQuad, screenparms.x0, 0);
screenparms.bottomLeftQuad  = OffsetRect(screenparms.topLeftQuad, 0, screenparms.y0);
screenparms.bottomRightQuad = OffsetRect(screenparms.topLeftQuad, screenparms.x0, screenparms.y0);


% ListenChar(2);  % Tell GetChar to start listening for keypresses
ShowCursor(0);	% arrow cursor
HideCursor;

screenparms.white = WhiteIndex(screenparms.window);
screenparms.black = BlackIndex(screenparms.window);

% Choose fonts likely to be installed on this platform
switch computer
    case 'MAC2',
        screenparms.serifFont     = 'Bookman';
        screenparms.sansSerifFont = 'Arial'; % or Helvetica
        screenparms.symbolFont    = 'Symbol';
        screenparms.displayFont   = 'Impact';
    case 'PCWIN'
        screenparms.serifFont     = 'Bookman Old Style';
        screenparms.sansSerifFont = 'Arial';
        screenparms.symbolFont    = 'Symbol';
        screenparms.displayFont   = 'Impact';
    case 'PCWIN64'
        screenparms.serifFont     = 'Bookman Old Style';
        screenparms.sansSerifFont = 'Arial';
        screenparms.symbolFont    = 'Symbol';
        screenparms.displayFont   = 'Impact';
    case 'GLNXA64'
        screenparms.serifFont     = 'Bookman Old Style';
        screenparms.sansSerifFont = 'Arial';
        screenparms.symbolFont    = 'Symbol';
        screenparms.displayFont   = 'Impact';
    otherwise
        error(['Unsupported OS: ' computer]);
end