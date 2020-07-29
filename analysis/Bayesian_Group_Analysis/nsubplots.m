function [nr, nc] = nsubplots(n)

% Determine ideal number of subplots for a given number of variables n
x{1} = round(sqrt(n));
x{2} = ceil(n/x{1});
[nr, nc] = deal(x{:});