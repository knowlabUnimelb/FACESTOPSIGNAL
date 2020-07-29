% Given a specific set of item conditions, extract the hits and false alarms
% This script will return two matrices:
%   The first is FA, false alarms rates from the same conditions (e.g., top
%   == 1).  
%   The second is the HR, hit rates from the different conditions (e.g.,
%   top ~= 1).
%   The columns are: top condition, bottom condition, response deadline,
%   false alarm rate or hit rate, trial count, total N

function [HR,FA,cols] = extractHandFA(si)

% optargs = {[1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]};
% newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
% optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
% [extract] = optargs{:}; % Place optional args in memorable variable names

%% Read in data
% si = 1; % which subject to analyse? [Must correspond to the entry in subjectNumber/subj variable
subjectNumbers = [2, 3, 4, 2, 3, 4]; % Subject number to analyse
sessions = {2:15, 2:15, 2:15, 2:15, 2:15, 2:15};
conditions = [1, 1, 1, 2, 2, 2];

subjectNumber = subjectNumbers(si);
session = sessions{si};
condition = conditions(si); 

useLongerDeadline = false; % Set to true to use 400 msec deadline

% Data columns
cols = {'subj', 'sess', 'blk', 'trial', 'respDuration',...
    'typeTop', 'typeBot', 'sameFlag', 'response', 'correctFlag', 'rt',...
    'rtFlag', 'stimCol1', 'stimCol2', 'stimCol3', 'stimCol4', 'stimCol5',...
    'stimCol6', 'stimCol7', 'stimCol8', 'stimCol9', 'stimCol10'};

% Read data
dataPrefix = '2019_facestopsignal';
[a,b] = strtok(fliplr(pwd), '\');
mainFolder = fliplr(b);
% datalocation = fullfile(mainFolder, '2019 FACESTOPSIGNAL', 'Data', 'Raw Data');
datalocation = fullfile(parentFolder(pwd, 2), 'data');
data = [];
for i = 1:numel(session)
    data = [data; dlmread(fullfile(datalocation, sprintf('%s_s%03d_con%d_ses%d.dat', dataPrefix, subjectNumber, condition, session(i))))];
end

%% Count good responses at each response deadline
newRTflag = data(:,strcmp(cols,'rtFlag')); 
newRTflag(isnan(newRTflag)) = 0;
pGood = aggregate([data(:,strcmp(cols, 'respDuration')), newRTflag], 1, 2, @nanmean);

%% Count response within 100 - 400 ms to see what happens if we increase the deadline
longRTflag = data(:,strcmp(cols,'rtFlag')); 
longRTflag(isnan(longRTflag)) = 0;
longRTflag(data(:,strcmp(cols, 'rt')) > 100 & data(:,strcmp(cols, 'rt')) < 400) = 1; 
pAlmostGood = aggregate([data(:,strcmp(cols, 'respDuration')), longRTflag], 1, 2, @nanmean);

%% Filter data based on rtFlag
pNans = mean(isnan(data(:,strcmp(cols, 'rtFlag'))));

longRTflag(isnan(data(:,strcmp(cols, 'rtFlag'))), :) = [];
data(isnan(data(:,strcmp(cols, 'rtFlag'))), :) = []; % This removes all of the nan RT flags, this is the response before onset of response signal

pGoodResponses = nanmean(data(:,strcmp(cols, 'rtFlag')));

switch useLongerDeadline
    case false
        data(data(:,strcmp(cols, 'rtFlag'))~=1, :) = []; % this removes all of the 0 RT flag (eg RT < 100 > 300) for all rows and columns
    case true
        data(~longRTflag,:) = [];
end

%% Get overall mean accuracy, standard deviations, counts, and convert to standard error
meanAcc = aggregate(data, mstrfind(cols, {'typeTop', 'typeBot', 'respDuration'}), mstrfind(cols, {'correctFlag', 'rt'}), @nanmean); 
stdsAcc = aggregate(data, mstrfind(cols, {'typeTop', 'typeBot', 'respDuration'}), mstrfind(cols, {'correctFlag', 'rt'}), @nanstd); 
cntsAcc = aggregate(data, mstrfind(cols, {'typeTop', 'typeBot', 'respDuration'}), mstrfind(cols, {'correctFlag', 'rt'}), @count); 
seAcc = stdsAcc(:,4)./sqrt(cntsAcc(:,4));

% Response deadlines are rounded oddly, fix it
meanAcc(:,3) = round(meanAcc(:,3)*100)./100;
stdsAcc(:,3) = round(stdsAcc(:,3)*100)./100;
cntsAcc(:,3) = round(cntsAcc(:,3)*100)./100;

%% Get HR and FA
FA = [];
FA = [meanAcc(meanAcc(:,1) == 1, 1:3), 1-meanAcc(meanAcc(:,1) == 1, 4), (1-meanAcc(meanAcc(:,1) == 1, 4)) .* cntsAcc(meanAcc(:,1) == 1, 4), cntsAcc(meanAcc(:,1) == 1, 4), meanAcc(meanAcc(:,1) == 1, 5)];
HR = [];
HR = [meanAcc(meanAcc(:,1) ~= 1, 1:3), meanAcc(meanAcc(:,1) ~= 1, 4), meanAcc(meanAcc(:,1) ~= 1, 4) .* cntsAcc(meanAcc(:,1) ~= 1, 4), cntsAcc(meanAcc(:,1) ~= 1, 4), meanAcc(meanAcc(:,1) ~= 1,5)];
cols = {'top', 'bot', 'dline', 'rate', 'Nresp', 'Ntot', 'rt'};
FA(:, strcmp(cols, 'Nresp')) = round(FA(:, strcmp(cols, 'Nresp')) );
HR(:, strcmp(cols, 'Nresp')) = round(HR(:, strcmp(cols, 'Nresp')) );