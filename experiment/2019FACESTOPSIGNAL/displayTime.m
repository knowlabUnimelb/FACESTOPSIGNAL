function [hours, mins, secs] = displayTime(totalTime, varargin)
% totalTime should be in seconds

optargs = {true};

% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);

% now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
optargs(newVals) = varargin(newVals);


% Place optional args in memorable variable names
[displayOn] = optargs{:};

%% 
hours = floor(totalTime/60/60);
% mins = floor(totalTime/60);
mins = floor((totalTime - hours * 60 * 60)/60);
% secs = round(totalTime - mins * 60);
secs = round((totalTime - hours * 60 * 60) - (mins * 60));

if displayOn
    disp(sprintf('Total time = %d seconds, or %d hours, %d minutes and %d seconds', round(totalTime), hours, mins, secs))
end