function [ydot,stepBound] = termChemoattractant(t,ENData,schemeData)

    % Data from the schemeData
    n = schemeData.n;
    membranePoints = schemeData.membranePoints{n-1};
    surfArea = schemeData.surfArea(n-1);
    grid     = schemeData.grid;
    Rc       = schemeData.Rc;

    % schemeData.surfArea is updated. volumePressure is assumed to point
    % outward normal to cell membrane.
    sigma_vol = -schemeData.volumeConst*(surfArea/schemeData.restingSurfArea-1);
    % sigma_ext is applied by the signal from LEGI-BEN scheme
    [sigma_ext] = applyChemoData(schemeData);
    % Pressure in cell
    % sigma_int = schemeData.pressureInCell;  %NOW PART OF sigma_ten
    % Model coefficients
    gamTen = schemeData.surfTension;
    km     = schemeData.cortexSpringCoefficient;
    gamm   = schemeData.cortexDampingCoefficient;
    gamc   = schemeData.anchorDampingCoefficient;

    % Input-output data
    phi = reshape(ENData(:,1), grid.shape);
    L   = reshape(ENData(:,2), grid.shape);

    % Compute curvature
    t1 = parameterizeByArcLength([membranePoints membranePoints membranePoints]);
    sp = csaps(t1,[membranePoints membranePoints membranePoints],0.98);
    t2 = linspace(min(t1),max(t1),3000);
    t2 = t2(1001:2001);
    x2 = fnval(sp,t2);
    dsp   = fnder(sp);   % First derivative?
    dspt  = fnval(dsp,t2); 
    ddspt = fnval(fnder(dsp),t2);  % Second derivative?
    % Curvature?    (x'y''-y'x'')/(x'^2+y'^2)^(3/2)
    kappa = -(dspt(1,:).*ddspt(2,:)-dspt(2,:).*ddspt(1,:))./(sum(dspt.^2)).^(3/2); %
    kappa = max(-100,min(kappa,100));
    kappa = griddata(x2(1,:),x2(2,:), kappa, grid.xs{1}, grid.xs{2}, 'nearest');
        % kappa = h.Feval('griddata',1,x2(1,:),x2(2,:), kappa, grid.xs{1}, grid.xs{2}, 'nearest');%flori
        % kappa = kappa{1};%flori

    % The cortical tension
    sigma_ten = gamTen.*(kappa-1/Rc);  % includes the old sigma_int

    % The total pressure
    sigma_tot = sigma_ext + sigma_vol - sigma_ten;
%     figure(1),surf(sigma_ext),shading flat, disp('sigma_ext'),[min(min(sigma_ext)),max(max(sigma_ext))]
%     disp('sigma_vol'),[min(min(sigma_vol)),max(max(sigma_vol))]
%     disp('sigma_int'),[min(min(sigma_int)),max(max(sigma_int))]
%     figure(2),surf(sigma_ten),shading flat, disp('sigma_ten'),[min(min(sigma_ten)),max(max(sigma_ten))]
%     figure(3),surf(sigma_tot),shading flat, 

     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normal velocity  
%     figure(4),surf(L),,shading flat,disp('L'),[min(min(L)),max(max(L))]
    speed = -(km./gamm) .* L + (1/gamc + 1./gamm) .* (sigma_tot);
%     figure(5),surf(speed),shading flat, disp('speed'),[min(min(speed)),max(max(speed))],
%          fprintf('\n sigma_ext %g',[max(max(sigma_ext))])
 
%     figure(6),surf(phi),shading flat, disp('phi'),[min(min(phi)),max(max(phi))],

    % All we care is the magnitude of grad(phi) and the stepBound
    magnitude    = zeros(grid.shape);
    stepBoundInv = zeros(grid.shape);

    for i = 1 : grid.dim
        % Get upwinded derivative approximations.
        [derivL,derivR] = feval(schemeData.derivFunc, grid, phi, i);
        % Effective velocity in this dimension (scaled by \|\grad \phi\|).
        prodL = speed .* derivL;
        prodR = speed .* derivR;
        magL = abs(prodL);
        magR = abs(prodR);
        % Determine the upwind direction.
        %   Either both sides agree in sign (take direction in which they agree),
        %   or characteristics are converging (take larger magnitude direction).
        flowL = ((prodL >= 0) & (prodR >= 0)) | ...
            ((prodL >= 0) & (prodR <= 0) & (magL >= magR));
        flowR = ((prodL <= 0) & (prodR <= 0)) | ...
            ((prodL >= 0) & (prodR <= 0) & (magL < magR));

        % For diverging characteristics, take gradient = 0
        %   (so we don't actually need to calculate this term).
        %flow0 = ((prodL <= 0) & (prodR >= 0));

        % Now we know the upwind direction, add its contribution to \|\grad \phi\|.
        magnitude = magnitude + derivL.^2 .* flowL + derivR.^2 .* flowR;

        % CFL condition: sum of effective velocities from O&F (6.2).
        effectiveVelocity = magL .* flowL + magR .* flowR;
        dxInv = 1 / grid.dx(i);
        stepBoundInv = stepBoundInv + dxInv * effectiveVelocity;
    end

    magnitude = sqrt(magnitude);
    clc
     fprintf('\n local radius %g %g',[min(abs(1./kappa(:))) max(abs(1./kappa(:)))])
     fprintf('\n Area %g',surfArea)
     fprintf('\n sigma_ext %g',[max(max(abs(sigma_ext)))])
     fprintf('\n sigma_vol %g',[max(max(abs(sigma_vol)))])
     fprintf('\n sigma_ten %g',[max(max(abs(sigma_ten)))])
     fprintf('\n sigma_tot %g',[max(max(abs(sigma_tot)))])
     fprintf('\n speed %g',[max(max(abs(speed)))])
     fprintf('\n deriv %g',max(max(max(abs(derivL))),max(max(abs(derivR)))))    
     fprintf('\n magnitude %g',[max(max(abs(magnitude)))])    
     fprintf('\n stepBoundInv %g',[max(max(abs(stepBoundInv)))])

    % Find the most restrictive timestep bound.
    nonZero = find(magnitude > 0);
    stepBoundInvNonZero = stepBoundInv(nonZero) ./ magnitude(nonZero);

    % Finally, compute the delta
    delta = speed .* magnitude;
    % Compute phiDot.
    phiDot = -delta(:);

    % Compute the stepBound
    stepBound = 1 / (max(stepBoundInvNonZero(:)));

    % Update the spring status

    LDot = -(km./gamm).*L+(1./gamm).*sigma_tot;

    % Compute the ydot
    ydot = [phiDot(:), LDot(:)];
end %termChemoattractant5

function sParameter = parameterizeByArcLength(dataPoints)
    % function [sParameter, dataPoints] = parameterizeByArcLength(dataPoints)
    % given an mxn array 'dataPoints' whose columns are [x y [data]]' at 
    % each of 'n' locations along a line, this function returns sParameter 
    % and dataPoints, where dataPoints is just like the input, except 
    % trimmed of any repeated values.  sParameter is a 1xn vector whose 
    % 'i'th element holds the arc length parameter of the 'i'th
    % dataPoint

    % fprintf('starting parameterize, dataPoints has size (%g, %g)\n',
    % size(dataPoints,1), size(dataPoints,2));

    pointsWithWrap = [dataPoints(:,end) dataPoints];
    % The x,y difference between subsequent points
    pointDifferences = dataPoints - pointsWithWrap(:,1:(end-1));
    segmentLengths = sqrt(pointDifferences(1,:).^2 + pointDifferences(2,:).^2);

    % Get rid of the point corresponding to position(1) - position(end) so that
    %   segmentLengths(1) corresponds to the length between points 1 and 2.
    segmentLengths(1) = [];
    % duplicatedPointIndexes = find(segmentLengths < 0.05);
    % segmentLengths(duplicatedPointIndexes) = [];
    % dataPoints(:,duplicatedPointIndexes) = [];

    % Computational trick -- it works out so that the "location" of the first
    % point is 0, and the "location" of the second point is its distance from
    % the first point
    segmentLengths = [0 segmentLengths];

    % OK, so I'm obsessed with this no for loop thing.  Element 'i' of
    % contourLocation is the sum of the elements 1 through 'i' of
    % segementLengths
    sParameter = tril(ones(size(segmentLengths,2),size(segmentLengths,2)))...
        *segmentLengths';
    sParameter = sParameter';
end

function [externalPressure] = applyChemoData(schemeData)

    %  This function transfer the response signal from the EN into protrusive
    %  stress
    n = schemeData.n;
    g = schemeData.grid;
    ntheta = size(schemeData.externalSignal,2);
    membranePoints = schemeData.membranePoints{n-1};
    centerX        = schemeData.centerX(n-1);
    centerY        = schemeData.centerY(n-1);

    % Calculate the angle of membrane points to the center within range of ...
    %   [-pi pi]
    membraneToCenterX = membranePoints(1,:) - centerX;
    membraneToCenterY = membranePoints(2,:) - centerY;
    angleToCenter = atan2(membraneToCenterY,membraneToCenterX);

    % for ntheta points
    % Transfer the external signal within the range of [-pi pi] (The ...
    %   original range is [0 2*pi])
    externalSignal = schemeData.externalSignal;
    the0=round(ntheta/2);
    externalSignal = [externalSignal(ntheta-the0+1:end) externalSignal(1:the0)];

    % Interpolate the signal value into all membrane points
    signalAngle = [-1:(2)/(ntheta-1):1]*pi;

    % To make sure signalAngle lies into range [-pi pi]
    signalAngle(end) = signalAngle(end) + 0.000001;  %PAI was 0.002

    membraneSignal = interp1(signalAngle,externalSignal,angleToCenter,'nearest');


    % Transfer the signal to force (PS: here, we don't have any rules now....
    %   We first try the proportional rule.
    membranePressure = 35 * membraneSignal; 

    % The following code is trying to make the contraction force larger
    membraneMask    = find(membraneSignal < 0);
    membranePressure(membraneMask) = membranePressure(membraneMask) * 1.5;

    externalPressure = griddata(membranePoints(1,:), membranePoints(2,:),...
        membranePressure,g.xs{1},g.xs{2},'nearest');

    % h = schemeData.h;%flori
    % externalPressure = h.Feval('griddata',1,membranePoints(1,:), membranePoints(2,:),...
    %     membranePressure,g.xs{1},g.xs{2},'nearest');%flori
    % externalPressure = externalPressure{1};%flori

end %applyChemoData