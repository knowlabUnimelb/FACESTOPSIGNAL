function saveStimImages(condition, imageName)

resizeProportion = .6;
screenparms = prepexp([0 0 0]);

%% Face Stimulus Set
if ismember(condition, [3 4])
    basebot = {'DX0101.bmp', 'DX0202.bmp', 'DX0303.bmp'};
    basetop = {'DY0101.bmp', 'DY0202.bmp', 'DY0303.bmp'};
elseif ismember(condition, [1 2])
    basetop = {'dDX0202.bmp', 'dDX0303.bmp', 'dDX0505.bmp'};
    basebot = {'dDY0202.bmp', 'dDY0303.bmp', 'dDY0505.bmp'};
else
    basetop = {'SX0101.bmp', 'SX0202.bmp', 'SX0303.bmp'};
    basebot = {'SY0101.bmp', 'SY0202.bmp', 'SY0303.bmp'};
end

% Define stimulus coordinates by the top half and the bottom half
stimuli = [3 3; 3 2; 2 3; 2 2; 3 1; 2 1; 1 3; 1 2; 1 1];

stimHeight = 400;
stimWidth  = 340;
baseRect   = [0 0 stimWidth stimHeight];
screenRect = [0 0 1280 1024];
stimRect   = CenterRect(baseRect, screenRect);

% These will vary by condition
% topRect    = [stimRect(1), stimRect(2), stimRect(3), stimRect(4)/2];
% botRect    = [stimRect(1), stimRect(4)/2, stimRect(3), stimRect(4)];

% Gaussian Blur vars
X = -stimWidth/2:stimWidth/2 - 1;
Y = -stimHeight/2:stimHeight/2 - 1;
XY = allcomb(X, Y);

if condition <= 4 % Morphs
    mX = 10; 
    mY = 1;
    sdX = 2000;
    sdY = 9000;
else              % Schematic Faces
    mX = 1;
    mY = 1;
    sdX = 3000;
    sdY = 9000;  
end
Z = mvnpdf(XY, [mX mY], [sdX, 0; 0 sdY]); % Gaussian Blur [Xmean, Ymean], [Xstd, 0; 0, Ystd]

% Bar separating top and bottom halves
barWidth   = 10;

% Offset
offsetWidth = 100;

stimTexture = cell(size(stimuli,1),1); % Initialize cell matrix for textures
for i = 1:size(stimuli,1)
    topInfo  = imfinfo(fullfile(pwd, 'BaseFaces', basetop{stimuli(i,1)}));
    botInfo  = imfinfo(fullfile(pwd, 'BaseFaces', basebot{stimuli(i,2)}));
    
    topImage = imread(fullfile(pwd, 'BaseFaces', basetop{stimuli(i,1)}), topInfo.Format);
    botImage = imread(fullfile(pwd, 'BaseFaces', basebot{stimuli(i,2)}), botInfo.Format);
    
    if size(topImage, 3) == 1; topImage = ind2rgb(imread(fullfile(pwd, 'BaseFaces',  basetop{stimuli(i,1)}), topInfo.Format), topInfo.Colormap) * 255; end
    if size(botImage, 3) == 1; botImage = ind2rgb(imread(fullfile(pwd, 'BaseFaces',  basebot{stimuli(i,2)}), botInfo.Format), botInfo.Colormap) * 255; end
    
    % If aligned then align top and bottom halves
    if ismember(condition, [1 3 5 7]);
        stimImage = zeros(size(topImage, 1), size(topImage, 2), 3);
        for j = 1:3
            stimImage(:,:,j) = [topImage(1:size(topImage,1)/2, :, j);
                                botImage(1 + size(topImage,1)/2:end, :, j)];
        end
    else
        stimImage = zeros(size(topImage, 1), offsetWidth + size(topImage, 2), 3);
        for j = 1:3
            stimImage(:,:,j) = [topImage(1:size(topImage,1)/2, :, j), zeros(size(topImage,1)/2, offsetWidth);
                                zeros(size(topImage,1)/2, offsetWidth), botImage(1 + size(topImage,1)/2:end, :, j)];
        end
    end
    
    % Add Gaussian blur
    Z3 =  repmat(reshape(Z, stimHeight, stimWidth), [1 1 3]); 
    if ~ismember(condition, [1 3 5 7])
        expandZ3 = nan(size(stimImage));
        for j = 1:3
            expandZ3(:, :, j) = [Z3(1:size(topImage,1)/2, :, j), zeros(size(topImage,1)/2, offsetWidth);
                                 zeros(size(topImage,1)/2, offsetWidth), Z3(1 + size(topImage,1)/2:end, :, j)];
        end
        Z3 = expandZ3;
    end
    stimImage = (stimImage .* Z3)./max(max(max((stimImage .* Z3)))) * 255;
    
    % Add black bar across split
    if condition <= 4
        stimImage(size(stimImage, 1)/2 - barWidth/2:size(stimImage, 1)/2 + barWidth/2, :, :) = 0;
    else
        stimImage(size(stimImage, 1)/2 - barWidth/2 + 10:size(stimImage, 1)/2 + barWidth/2 + 10, :, :) = 0;
    end
    
    % Resize 
    stimImage = imresize(stimImage, resizeProportion, 'Colormap', 'original');
    if resizeProportion < 1
        stimImage(stimImage < 0 ) = 0; stimImage(stimImage > 255) = 255; 
    end
    
    
    % Get rect
    [stimRect,dh,dv] = CenterRect([0 0 size(stimImage, 1), size(stimImage, 2)], screenparms.rect);
    
    % Save to texture
    stimTexture{i} = Screen('MakeTexture', screenparms.window, stimImage);
end



for i = 1:numel(stimTexture)
    if ismember(condition, [1 2 5 6]);
        Screen('DrawTexture', screenparms.window, stimTexture{i});
    else
        Screen('DrawTexture', screenparms.window, stimTexture{i}, [], [], 180);
    end
    image = Screen('GetImage', screenparms.window, screenparms.rect, 'backBuffer');
    imwrite(image, fullfile(pwd, sprintf('%s_%d.bmp', imageName, i)));
end
closeexp(screenparms);