function schemeData = Initialization(schemeData, membranePotentialFunction,levelSetConfiguration)

%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%

g = schemeData.grid;
Rc = schemeData.Rc;

%Find the 0 level set as the initial membrane boundary
%   (c = contourc( x, y, z, v )) where v is a vector
membranePoints = contourc(g.xs{1}(:,1),g.xs{2}(1,:),...
    membranePotentialFunction', [0 0]);

%Remove the first column, since it only contains level and number of pairs  
membranePoints = membranePoints(:,2:(membranePoints(2,1)+1));

if (polygon_area_2d(membranePoints)>=0)  %if curve is ccw, make it clockwise
    membranePoints = membranePoints(:,end:-1:1);
end

% The stats for cell
[geom, iner, cpmo] = polygeom(membranePoints(1,:)',membranePoints(2,:)');

% Record the spring state
sizeMembrane = size(membranePoints);
numMemPoints = sizeMembrane(2);
memSpring    = zeros(1,numMemPoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial schemeData parameters (ideal and resting circumference of circle) 
%  for future plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

schemeData.membranePoints{1} = membranePoints;
schemeData.memSpring{1} = memSpring;
schemeData.surfArea(1)  = geom(1);
schemeData.perimeter(1) = geom(4);
schemeData.centerX(1)   = geom(2);
schemeData.centerY(1)   = geom(3);
schemeData.time(1)      = 0;
schemeData.n            = 2;
schemeData.levelSetConfiguration = levelSetConfiguration;