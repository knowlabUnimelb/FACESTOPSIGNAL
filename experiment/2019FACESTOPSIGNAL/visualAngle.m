function angle = visualAngle(viewingdistance, stimulussize)
% VISUALANGLE(D, S) Compute the visual angle subtended by a stimulus of size S
% given a viewing distance of D
%   
% distances should be in centimeters (though it doesn't matter as long as
% they're consistent)

%%  Variables
adj = viewingdistance;
opp = stimulussize;

%% Dan's formula
% hyp = sqrt(adj^2 + opp^2); % Hypotenuse
% 
% 
%angle = degrees(asin(opp/hyp)) * 2; % Angle

%% Wiki's formula

angle = 2*radtodeg(atan(opp/(2*adj)));

end
