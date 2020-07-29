function output = parentFolder(folder, varargin)

optargs = {1};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[nMoves] = optargs{:}; % Place optional args in memorable variable names

output = fileparts(folder);
idx = 1;
while idx < nMoves
    output = fileparts(output);
    idx = idx + 1; 
end    
