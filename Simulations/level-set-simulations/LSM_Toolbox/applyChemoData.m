function externalPressure = applyChemoData(schemeData)

% function externalPressure = applyChemoData(schemeData);
%
%  This function transfer the response signal into pressure

%%%%%%%%%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ?schemeData
%%%%%%%%%Output%%%%%%%%%%%%%
% externalPressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the transformation law is not clear now, we just make it up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = schemeData.n;
g = schemeData.grid;
membranePoints = schemeData.membranePoints{n-1};
centerX = schemeData.centerX(n-1);
centerY = schemeData.centerY(n-1);

% Calculate the angle of membrane points to the center within range of ...
%   [-pi pi]
membraneToCenterX = membranePoints(1,:) - centerX;
membraneToCenterY = membranePoints(2,:) - centerY;
angleToCenter = atan2(membraneToCenterY,membraneToCenterX);

% % for 720 points
% % Transfer the external signal within the range of [-pi pi] (The ...
% %   original range is [0 2*pi])
% externalSignal = schemeData.externalSignal;
% externalSignal = [externalSignal(362:end) externalSignal(1:361)];
% 
% % Interpolate the signal value into all membrane points
% signalAngle = [-pi:(2*pi)/719:pi]; 

% for 314 points
% Transfer the external signal within the range of [-pi pi] (The ...
%   original range is [0 2*pi])
externalSignal = schemeData.externalSignal;
externalSignal = [externalSignal(159:end) externalSignal(1:158)];

% Interpolate the signal value into all membrane points
signalAngle = [-1:(2)/313:1]*pi; 

% To make sure signalAngle lies into range [-pi pi]
signalAngle(end) = signalAngle(end) + 0.002;

membraneSignal = interp1(signalAngle,externalSignal,angleToCenter,'nearest');

% Transfer the signal to force (PS: here, we don't have any rules now....
%   We first try the proportional rule.
membranePressure = 35 * membraneSignal; % Force should be in the rang of 1-2 pN/um2 max

% The following codes are trying to make the contraction force larger
% membraneMask = find(membraneSignal < 0);
% membranePressure(membraneMask) = membranePressure(membraneMask) * 1.5;

externalPressure = griddata(membranePoints(1,:), membranePoints(2,:),...
    membranePressure,g.xs{1},g.xs{2},'nearest');
% h = schemeData.h;%flori
% externalPressure = h.Feval('griddata',1,membranePoints(1,:), membranePoints(2,:),...
%     membranePressure,g.xs{1},g.xs{2},'nearest');%flori
% externalPressure = externalPressure{1};%flori
