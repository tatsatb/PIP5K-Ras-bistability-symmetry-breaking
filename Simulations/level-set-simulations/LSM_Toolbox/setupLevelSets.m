function levelSetConfiguration = setupNewLevelSets(accuracy,range,gs,volConst,surfTen,visco,Rc);

%function levelSetConfiguration = setupNewLevelSets(accuracy,range,gs,volConst,surfTen,visco,Rc);
%
% This function creates a variable levelSetConfiguration which holds all of 
% the constants and settings required to do a level set simulation.
%
% Inputs  -----------------------------------------------------------------
%
% yaccuracy:    low, medium, high, veryHigh - how good of a solver do you
%               want to use to solve your level sets
% range:        [xleft, xright, ybottom, ytop]
%               xleft, xright - the range of x values over which you wish 
%               to create a simulation grid
%               ybottom, ytop - the range of y values over which you wish 
%               to create a simulation grid
% gs:           grid spacing 
% volConst:     Volume Conservation Constant;
% surfTen:      Surface Tension constant
% visco:        Viscoelastic parameters (K, D, B as in Yang, 2008)
% Rc:           Initial radius of cell (in um)
%
% Outputs  ----------------------------------------------------------------
%
% levelSetConfiguration: variable holding all of the configuration data.
%   .integratorFunc:    which of the integrators included in the level
%                       set toolbox should be use to integrate and solve 
%                       the differential equations
%   .integratorOptions: Any options required for integratorFunc
%   .schemeData
%       .grid:  the grid structure that defines the area over which the
%               simulation occurs (both in terms of bounds and resolution).
%       .xs{1}- all the y-values of the grid in columns
%       .xs(2)- all the x-values of the grid in rows
%       .areaB- The gain constant in the P controller that enforces
%               the conservation of area.
%       .anchorMembraneSpringCoefficient: scalar that determines how 
%               much force is exerted by the membrane if it is stretched 
%               (e.g. Hooke's law: F = k*x).
%       .anchorDampingCoefficient
%       .membraneDampingCoefficient

xleft   = range(1);
xright  = range(2); 
ybottom = range(3);  
ytop    = range(4);
K       = visco(1);
D       = visco(2);
B       = visco(3);

if(nargin ~= 1)
    accuracy = 'low';
end

% Create the grid.
g.dim = 2;
g.min = [xleft;  ybottom];
g.max = [xright; ytop];
g.N = [((xright-xleft)*gs+1); ((ytop-ybottom)*gs+1)];
g.bdry = @addGhostExtrapolate;  % A ToolBox function

% This fills in all the missing parts of grid.___ structure based on the above dimensions
g = processGrid(g);
schemeData.grid = g;
schemeData.volumeConst = volConst;
schemeData.surfTension = surfTen;
schemeData.cortexSpringCoefficient  = K;
schemeData.cortexDampingCoefficient = D; 
schemeData.anchorDampingCoefficient = B;
schemeData.Rc = Rc;
schemeData.restingSurfArea = pi*Rc^2;
schemeData.restingPerimeter = 2*pi*Rc;
schemeData.pressureInCell = surfTen*2/Rc;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'off');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
case 'low'
 schemeData.derivFunc = @upwindFirstFirst;
 integratorFunc = @odeCFL1Shi;
case 'medium'
 schemeData.derivFunc = @upwindFirstENO2;
 integratorFunc = @odeCFL2;
case 'high'
 schemeData.derivFunc = @upwindFirstENO3;
 integratorFunc = @odeCFL3;
case 'veryHigh'
 schemeData.derivFunc = @upwindFirstWENO5;
 integratorFunc = @odeCFL3;
otherwise
 error('Unknown accuracy level %s', accuracy);
end

levelSetConfiguration.schemeData = schemeData;
levelSetConfiguration.integratorFunc = integratorFunc;
levelSetConfiguration.integratorOptions = integratorOptions;

