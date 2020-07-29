function showInstructions(screenparms, instructionFile, advance, waittime)
global ptb3
ptb3 = true;

if nargin == 2
    advance = 'space';
end

imageInfo = imfinfo(instructionFile);
instructionImage = imread(instructionFile, imageInfo.Format);
if size(instructionImage, 3) == 1
    instructionImage = ind2rgb(imread(instructionFile, imageInfo.Format), imageInfo.Colormap) * 255;
end


[imageRect,dh,dv] = CenterRect([0 0 imageInfo.Width imageInfo.Height], screenparms.rect);
if screenparms.rect(3)/imageInfo.Width < 1 || screenparms.rect(4)/imageInfo.Height < 1
    newWidth = imageInfo.Width * min([screenparms.rect(3)/imageInfo.Width, screenparms.rect(4)/imageInfo.Height]);
    newHeight = imageInfo.Height * min([screenparms.rect(3)/imageInfo.Width, screenparms.rect(4)/imageInfo.Height]);
    [imageRect,dh,dv] = CenterRect([0 0 newWidth newHeight], screenparms.rect);
end


if ptb3
    instructionTexture = Screen('MakeTexture', screenparms.window, instructionImage);
    Screen('DrawTexture', screenparms.window, instructionTexture, [], imageRect);
    Screen('Flip', screenparms.window);
else
    Screen(screenparms.window, 'PutImage',  instructionImage, imageRect);
end

switch advance
    case 'space'
        try
            pressSpace; if ptb3; Screen('Flip', screenparms.window); end
        catch
            pause; FillScreen(screenparms);  if ptb3; Screen('Flip', screenparms.window); end
        end
    case 'time'
        if nargin == 3
            waittime = 1;
        end
        WaitSecs(waittime);
    case 'RTBox'
          RTBox('clear'); % clear buffer and sync clocks before stimulus onset
          while ~any(RTBox('ButtonDown')); WaitSecs(0.01); end; 
          WaitSecs(.1);
end
if ptb3; Screen('Close', instructionTexture); end 