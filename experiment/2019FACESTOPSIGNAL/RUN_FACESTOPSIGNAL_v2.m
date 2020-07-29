%% Face Stop signal task
clear variables
clc

WaitSecs(1e-7); % Hack to load WaitSecs
warning('off', 'MATLAB:mex:deprecatedExtension')

debug = false; % If debug == true, do experiment with reduced number of trials
RTBox('fake', debug)

debug = false; % will set transparent screen for debugging
if debug % The statement below will run if the value of study.debug_mode = 1
    PsychDebugWindowConfiguration(0.5); % This produces a transparent screen so we can see what's happening on the command line while the experiment is running
    % This skips a whole bunch of tests which are normally done when we use a CRT monitor or are running the real experiment
end

%% Experimental Variables
seed = 41283;
xup = -50; xleft = -18; % center of screen offsets for showtext
bgcolor = [0 0 0]; % User specified background color
txtColor = [255 255 255];

textsize = 30;

fixsize   = 50;  % Height and width of fixation cross (in pixels)
lineWidth = 5;   % Pixel width

resizeProportion = .6;

fixationDuration = .5;  % Fixation cross presentation length (1.5 secs)
postProbeDuration = .5;
fbkAccDuration   = .25; % Accuracy Feedback presentation length
postFbkAccBlank  = .4;  % Post accuracy feedback blank interval
fbkRTDuration    = .5; 
postFbkRTBlank   = .3; 
fbkTooEarly = 2; 
timeout = 5;

%% Trial numbers
nExpTrials = 780; % Number of experimental trials = 10 blocks of 78 trials
nBlocks    = 10;  % Number of blocks = 10
nTrialsPerBlock = nExpTrials/nBlocks; % Number of trials per block = 78

responseDurations = [.050, .100, .200, .400, .800, 1.800];

probeDuration = .5; 
rtTimeOut = .3; % RT must be within 300 msec after response signal
rtTooFast = .1; % RT must be after the response signal

%% Open Experiment Window
% Present subject information screen until correct information is entered
subject     = input('Enter Subject Number [101-110]:');
session     = input('Enter Session Number [2-15]:');
condition   = input('Enter Condition Number [1 or 2]:'); % 1 for aligned, 2 for misaligned

% screenparms = prepexp(0, bgcolor, [0 0 1280 1024]); % Open onscreen window
screenparms = prepexp([], bgcolor); % Open onscreen window
rng(seed + session * subject * 2)

fprintf('Window loaded\n')
outputfile = sprintf('2019_facestopsignal_s%03d_con%d_ses%d.dat', subject, condition, session);

%% Face Stimulus Set
facesets = 1:20; % Each trial will need to grab a top and bottom from a different face set (probe and test faces, tops and bottoms, will be drawn from the same sets)

% Create large vector of possible test images for each base image
% Test images for the top will be: same, same, same, low, high, low, low, high, high
% Test images for the bot will be: same, low, high, same, same, low, high, low, high
testsets = [1 1; 1 2; 1 3; 2 1; 3 1; 2 2; 2 3; 3 2; 3 3];

% Probabilities for different trial types
testFaceProbabilities = [1/6 1/6 1/6 1/12 1/12 1/12 1/12 1/12 1/12];

% Same, L = Hard, H = Easy discriminations are based on the base face, each row indicates an increment referenced in test sets
testidx  = {[1 1; 2 2; 3 3], [1 3; 2 2; 3 1], [], [3 1; 2 2; 1 3], [3 3; 2 2; 1 1]}; 

% These are the indexes which specify how the images are loaded
allidx   = [1 1; 1 3; 2 2; 3 1; 3 3]; % These are references to individual faces
allStimIdx = allcomb(facesets, 1:size(allidx, 1));
fullStimIdx = [allStimIdx(:,1), allidx(allStimIdx(:,2), :)]; 
%%
stimidx = allStimIdx(allStimIdx(:,2) ~= 3, :);

%% LOAD STIMULI INTO TEXTURES
% The easiest/least memory requirements way to do this is to load each
% face, then present them on screen in the appropriate position with the
% bar and blurring applied after the fact; I'll try blurring the faces first
% can I present half an image?

% Initialize progress bar
ProgBar = DrawProgressBar(screenparms, size(allStimIdx,1), 'Generating Stimuli', txtColor); Screen('Flip', screenparms.window);

stimHeight = 400;
stimWidth  = 340;
baseRect   = [0 0 stimWidth stimHeight];
screenRect = [0 0 1280 1024];
% stimRect   = CenterRect(baseRect, screenRect);

% Gaussian Blur vars
X = -stimWidth/2:stimWidth/2 - 1;
Y = -stimHeight/2:stimHeight/2 - 1;
XY = allcomb(X, Y);

mX = 2;
mY = -10;
sdX = 2000;
sdY = 8000;

Z = mvnpdf(XY, [mX mY], [sdX, 0; 0 sdY]); % Gaussian Blur [Xmean, Ymean], [Xstd, 0; 0, Ystd]

% Bar separating top and bottom halves
barWidth   = 10;

% Offset
% offsetWidth = 00;
if condition == 1
    offsetWidth = 00;
else
    offsetWidth = 85; %added this in [Jun]
end

stimTexture = cell(size(allStimIdx,1),1); % Initialize cell matrix for textures
for i = 1:size(allStimIdx,1)
    faceName  = sprintf('ss_set%02d%02d%02d.bmp', facesets(allStimIdx(i,1)), allidx(allStimIdx(i,2), 1), allidx(allStimIdx(i,2), 2));
    faceInfo  = imfinfo(fullfile(pwd, 'FaceSets', faceName));
    %     faceImage = imread(fullfile(pwd, 'FaceSets', faceName), faceInfo.Format);
    faceImage = ind2rgb(imread(fullfile(pwd, 'FaceSets', faceName), faceInfo.Format), faceInfo.Colormap) * 255; % Need this to set the colourmap appropriately
    blurImage = faceImage;
    
    
    if condition == 1 % start of code I'm experimenting with [Jun]
        stimImage = zeros(size(blurImage, 1), size(blurImage, 2), 3);
        for j = 1:3
            stimImage(:,:,j) = [blurImage(1:size(blurImage,1)/2, :, j);
                                blurImage(1 + size(blurImage,1)/2:end, :, j)];
        end
    else
        stimImage = zeros(size(blurImage, 1), offsetWidth + size(blurImage, 2), 3);
        for j = 1:3
                stimImage(:,:,j) = [blurImage(1:size(blurImage,1)/2, :, j), zeros(size(blurImage,1)/2, offsetWidth);
                    zeros(size(blurImage,1)/2, offsetWidth), blurImage(1 + size(blurImage,1)/2:end, :, j)];
        end
    end
    
    % Add Gaussian blur [Dan's code]
    Z3 =  repmat(reshape(Z, stimHeight, stimWidth), [1 1 3]);
        
    expandZ3 = nan(size(blurImage, 1), offsetWidth + size(blurImage, 2), 3);
    for j = 1:3
        expandZ3(:, :, j) = [Z3(1:size(blurImage,1)/2, :, j), zeros(size(blurImage,1)/2, offsetWidth);
            zeros(size(blurImage,1)/2, offsetWidth), Z3(1 + size(blurImage,1)/2:end, :, j)];
    end
    Z3 = expandZ3;
    stimImage = (stimImage .* Z3)./max(max(max((stimImage .* Z3)))) * 255;

    % Add black bar across split
    stimImage(size(stimImage, 1)/2 - barWidth/2:size(stimImage, 1)/2 + barWidth/2, :, :) = 0; %end of code I'm experimenting with [Jun]
    
    % Resize
    blur2 = stimImage;
    stimImage = imresize(blur2, resizeProportion, 'Colormap', 'original');
    if resizeProportion < 1
        stimImage(stimImage < 0 ) = 0; stimImage(stimImage > 255) = 255;
    end
    
    % Save to texture
    stimTexture{i} = Screen('MakeTexture', screenparms.window, stimImage);
    ProgBar(i); Screen('Flip', screenparms.window); % Update the progress bar
end

% Get rect locations
faceRect = [0 0 size(stimImage, 2), size(stimImage, 1)];
[stimRect,dh,dv] = CenterRect(faceRect, screenparms.rect);
[w, h] = RectSize(stimRect);
topRectSource = [faceRect(1), faceRect(2), faceRect(3), faceRect(2)+h/2];
botRectSource = [faceRect(1), faceRect(2)+h/2, faceRect(3), faceRect(4)];
topRectDest = [stimRect(1), stimRect(2), stimRect(3), stimRect(2)+h/2];
botRectDest = [stimRect(1), stimRect(2)+h/2, stimRect(3), stimRect(4)];

% Set bar locations
barloc = CenterRect([0 0 size(stimImage, 2) barWidth], screenparms.rect);

% Fixation cross locations
fixcrossH = CenterRect([0 0 fixsize 1], screenparms.rect);
fixcrossV = CenterRect([0 0 1 fixsize], screenparms.rect);

%% PRESENT INSTRUCTIONS
instructionFolder = 'Instructions v2';
trainingInstructions = {'Instructions_Experiment.bmp'};
showInstructions(screenparms, fullfile(pwd, instructionFolder, trainingInstructions{1}), 'RTBox')
endImage = (fullfile(pwd, instructionFolder, 'Thanks.bmp'));

%% Load crashed data file?
if exist(outputfile, 'file') == 2
    tempdata =  dlmread(outputfile);
    lastCompletedBlock = tempdata(end,3);
    startBlock = lastCompletedBlock+1;
else
    startBlock = 1;
end


%% Run Experiment
attendToTop = 1; 
for bblocks = startBlock:nBlocks
    output = [];
    
    % Set up response durations
    rds = repmat(responseDurations, nTrialsPerBlock./numel(responseDurations), 1);
    responseDuration = rds(:);
    responseDuration = responseDuration(randperm(numel(responseDuration)));
    
    % Preallocate
    stimuli = nan(nTrialsPerBlock, 10);
    type = nan(nTrialsPerBlock, 2);
    rt = nan(nTrialsPerBlock, 1);
    response = nan(nTrialsPerBlock, 1);
    sameFlag = nan(nTrialsPerBlock, 1);
    correctFlag = nan(nTrialsPerBlock, 1);
    rtFlag = nan(nTrialsPerBlock, 1);
    
    % Start experiment
    priorityLevel = MaxPriority(screenparms.window,'WaitBlanking');
    trialcnt = 1;
    tic;
    for i = 1:nTrialsPerBlock

        % Select base faces
        top = stimidx(datasample(1:size(stimidx, 1), 1), :);
        bot = stimidx(datasample(1:size(stimidx, 1), 1), :);
        
        probeTopTex = find(all(allStimIdx == ones(size(allStimIdx, 1),1) * top, 2));
        probeBotTex = find(all(allStimIdx == ones(size(allStimIdx, 1),1) * bot, 2));
        
        % Select test face
        testid = datasample(1:numel(testFaceProbabilities), 1, 'Weights', testFaceProbabilities); % Determine what type of test trial
        testtype = testsets(testid,:); % Determines the change in the top and bottom: 1 = same, 2 = low, 3 = high
        type(i,:) = testtype;
        
        testTopTex = find(all(fullStimIdx == ones(size(fullStimIdx, 1),1) * [top(1) testidx{top(2)}(testtype(1), :)], 2));
        testBotTex = find(all(fullStimIdx == ones(size(fullStimIdx, 1),1) * [bot(1) testidx{bot(2)}(testtype(2), :)], 2));
        
        stimuli(i,:) = [fullStimIdx(allStimIdx(:,1) == top(1) & allStimIdx(:,2) == top(2), :),...
            testidx{top(2)}(testtype(1), :),...
            fullStimIdx(allStimIdx(:,1) == bot(1) & allStimIdx(:,2) == bot(2), :),...
            testidx{bot(2)}(testtype(2), :)];
        
        % Display fixation cross
        Screen('DrawLine', screenparms.window, screenparms.white, fixcrossH(1), fixcrossH(2), fixcrossH(3), fixcrossH(4), lineWidth)
        Screen('DrawLine', screenparms.window, screenparms.white, fixcrossV(1), fixcrossV(2), fixcrossV(3), fixcrossV(4), lineWidth)
        Screen('Flip', screenparms.window); % Flip the fixation cross to the front buffer
        WaitSecs(fixationDuration); % Wait
        
        Priority(priorityLevel);
        RTBox('clear', 10)
        
        % Present probe top and bottom of textures
        Screen('DrawTexture', screenparms.window, stimTexture{probeTopTex}, topRectSource, topRectDest);
        Screen('DrawTexture', screenparms.window, stimTexture{probeBotTex}, botRectSource, botRectDest);
        
        % Draw black bar
        Screen('FillRect', screenparms.window, [0 0 0], barloc)
        Screen('Flip', screenparms.window);
        WaitSecs(probeDuration);
        
        % Black screen
        FillScreen(screenparms);
        Screen('Flip', screenparms.window);
        WaitSecs(postProbeDuration);
        
        % Present test top and bottom of textures
        Screen('DrawTexture', screenparms.window, stimTexture{testTopTex}, topRectSource, topRectDest);
        Screen('DrawTexture', screenparms.window, stimTexture{testBotTex}, botRectSource, botRectDest);
        
        % Draw black bar
        Screen('FillRect', screenparms.window, [0 0 0], barloc)
        Screen('Flip', screenparms.window);
        
        % If a button is pressed before the onset of the response signal
        % stop in error and display penalty text
        % WaitSecs(responseDuration(i))
        startTime = GetSecs;
        earlyResponse = false;
        while (GetSecs - startTime) <= responseDuration(i)
            isDown = RTBox('ButtonDown');
            if any(isDown)
                % Stop trial and yell at participant for not waiting
                earlyResponse = true;
                break
            end
        end
        
        if ~earlyResponse
            % Present test face with response instructions
            Screen('DrawTexture', screenparms.window, stimTexture{testTopTex}, topRectSource, topRectDest);
            Screen('DrawTexture', screenparms.window, stimTexture{testBotTex}, botRectSource, botRectDest);
            Screen('FillRect', screenparms.window, [0 0 0], barloc)
            showtext(screenparms, textsize, 'xxxxxxxxxxxxxxxxxxxxxxxxx', 0, -h, -w/8);
            showtext(screenparms, textsize, 'xxxxxxxxxxxxxxxxxxxxxxxxx', 0,  h, -w/8);
            %showtext(screenparms, textsize, 'xxxxx RESPOND NOW xxxxx', 0, -h, -w/3);
            %showtext(screenparms, textsize, 'xxxxx RESPOND NOW xxxxx', 0,  h, -w/3);
            RTBox('clear'); % clear buffer and sync clocks before stimulus onset
            vbl = Screen('Flip', screenparms.window);
            
            
            % Record response
            [cpuTime, buttonPress] = RTBox(timeout);  % computer time of button response
            FillScreen(screenparms);
            Screen('Flip', screenparms.window);
            
            if ~isempty(cpuTime); rt(i, 1) = (cpuTime(1) - vbl) * 1000;
            else 
                rt(i, 1) = nan;
                buttonPress = '9';
            end
            
            if     ismember(buttonPress(1), {'1', '2'}); response(i, 1) = 1;   % Same
            elseif ismember(buttonPress(1), {'3', '4'}); response(i, 1) = 2;   % Diff
            else                                         response(i, 1) = nan;
            end
            
            % check response and rt
            switch attendToTop
                case 1
                    if testtype(1) == 1; sameFlag(i,1) = 1; % Top is the same
                    else                 sameFlag(i,1) = 2; % Top is different
                    end
                case 0
                    if testtype(2) == 1; sameFlag(i,1) = 1; % Bot is the same
                    else                 sameFlag(i,1) = 2; % Bot is different
                    end
            end
            
            if response(i,1) == sameFlag(i,1)
                correctFlag(i,1) = 1;
            else
                correctFlag(i,1) = 0;
            end
            
            if rt(i,1) < rtTimeOut*1000 && rt(i,1) > rtTooFast*1000
                rtFlag(i,1) = 1;
            elseif rt(i,1) < rtTimeOut*1000 && rt(i,1) < rtTooFast*1000
                rtFlag(i,1) = 2;
            else
                rtFlag(i,1) = 0;
            end
            
            % Display feedback
            % Feedback on accuracy for 250 ms, 400 ms blank
            if correctFlag(i,1)
                showtext(screenparms, textsize, 'CORRECT', 0, 0, -20);
                %showtext(screenparms, textsize, 'Response is Correct', 0, 0, 0);
            else
                showtext(screenparms, textsize, 'WRONG', 0, 0, -20);
                %showtext(screenparms, textsize, 'Response is Wrong', 0, 0, 0);
            end
            Screen('Flip', screenparms.window);
            WaitSecs(fbkAccDuration)
            
            Screen('Flip', screenparms.window);
            WaitSecs(postFbkAccBlank)
            
            if rtFlag(i,1) == 1
                showtext(screenparms, textsize, 'ON TIME', 0, 0, -20);
                %showtext(screenparms, textsize, 'Response Time: Good', 0, 0, 0);
            elseif rtFlag(i,1) == 2
                showtext(screenparms, textsize, 'TOO FAST', 0, 0, -20);
                %showtext(screenparms, textsize, 'Response Time: TOO FAST', 0, 0, 0);
            else
                showtext(screenparms, textsize, 'TOO SLOW', 0, 0, -20);
                %showtext(screenparms, textsize, 'Response Time: TOO SLOW', 0, 0, 0);
            end
            Screen('Flip', screenparms.window);
            WaitSecs(fbkRTDuration)
            
            Screen('Flip', screenparms.window);
            WaitSecs(postFbkRTBlank)
            Priority(0);
        else % Show responded too early feedback
            showtext(screenparms, textsize, 'PLEASE WAIT FOR THE RESPONSE SIGNAL', 0, 0, -20);
            Screen('Flip', screenparms.window);
            WaitSecs(fbkTooEarly)
        end
        
        %% Insert a short break at 1/4, 1/2, and 3/4 trials
        trialcnt = trialcnt + 1;
    end
  
    % Save output
    output = [output; 
              bblocks * ones(nTrialsPerBlock, 1),...
              (1:nTrialsPerBlock)',...
              responseDuration,...
              type,...
              sameFlag,...
              response,...
              correctFlag,...
              rt,...
              rtFlag,...
              stimuli]; 
    
    if bblocks == nBlocks
        showtext(screenparms, 20, 'Preparing Output', 0, 0, 0);
        Screen('Flip', screenparms.window); 
    else % Insert a break
        showtext(screenparms, 20, 'Take a short break. Press any button to continue', 0, 0, 0); % Show some instructions to advance after the break
        Screen('Flip', screenparms.window);
        
        RTBox('clear');                                                % Clear buffer and sync clocks before stimulus onset
        while ~any(RTBox('ButtonDown')); WaitSecs(0.01); end           % Wait for any button press
        Screen('Flip', screenparms.window);                            % Clear the screen after a button press        
    end
        
    
    %% Save data
    finaloutput = [subject * ones(trialcnt-1,1), session * ones(trialcnt-1,1), output];
    dlmwrite(outputfile, finaloutput, '-append');
    displayTime(toc)
end
showInstructions(screenparms, endImage, 'RTBox');
closeexp(screenparms) %% Close experiment

% Columns for reference
cols = {'subj', 'sess', 'blk', 'trial', 'respDuration', 'typeTop', 'typeBot', 'sameFlag', 'response', 'correctFlag', 'rt',    'rtFlag', 'stimCol1', 'stimCol2', 'stimCol3', 'stimCol4', 'stimCol5', 'stimCol6', 'stimCol7', 'stimCol8', 'stimCol9', 'stimCol10'};