function [charx, t] = getresponse

FlushEvents('keyDown');
startSecs = GetSecs;
while 1
    [keyIsDown,timeSecs,keyCode] = KbCheck;
    if keyCode(KbName('f12'))
        Screen('closeall'); FlushEvents('keyDown');
    end
    if keyIsDown
        c = KbName(keyCode);
        if iscell(c)
            charx = c{1}(1);
        else
            charx = c(1);
        end
        t = timeSecs - startSecs;
        while KbCheck; end
        break;
    end
end

FlushEvents('keyDown');