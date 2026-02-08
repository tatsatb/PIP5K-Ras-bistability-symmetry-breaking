function LSM_moving_cell(param,file_type,ext_Signal_norm,folder_name)
% Script to generate level set simulations using LSM toolbox

    pt = pwd;
    addPathToKernel(pt);




    %% folder to save the LS simulation output

    
    if ~exist(folder_name,'dir')
        mkdir(folder_name)
    end

    %%---------------------------------------------------------------------
    % Set up model kinetics and coefficients

    Rc  = param.Rc;   % Initial cell radius in um
    % The cell's mechanical properties; from Yang et al., BMC Sys Bio, 2008
    km   = param.km; % s*nN/um^2
    gamm = param.gamm; % s*nN/um^3
    gamc = param.gamc;  % s*nN/um^3
    visco= [km, gamm, gamc];

    volCKp = param.volCKp; % (dimensionless);
    volCKi = param.volCKi;  % (dimensionless)
    gamTen = param.gamTen; % nN/um; from Yang et al., BMC Sys Bio, 2008

    gs = param.gs; % Grid numbers per um 

    %----------------------------------------------------------------------
    % Set up level sets

    % Sets up grid dimensions automatically based on cell size.
    xleft_grid   = param.xleft_grid;
    xright_grid  = param.xright_grid;
    ybottom_grid = param.ybottom_grid;
    ytop_grid    = param.ytop_grid;
    palette      = [xleft_grid,xright_grid,ybottom_grid,ytop_grid];

    % SetupLevelSets returns a levelSetConfiguration which contains setup
    % parameters for simulation

    levSetConfiguration = setupLevelSets('low',palette,gs,volCKp,gamTen,visco,Rc);

    schemeData         = levSetConfiguration.schemeData;
    integratorFunc     = levSetConfiguration.integratorFunc;
    integratorOptions  = levSetConfiguration.integratorOptions;
    g = schemeData.grid;

    %----------------------------------------------------------------------
    % Initialization

    % fprintf('Creating Initial Potential Function ... ');
    membranePotentialFunction = shapeSphere(g, [0 0],Rc);
    % fprintf('done!\n')

    % fprintf('Initializing schemeData values ... ');
    schemeData = Initialization(schemeData, membranePotentialFunction,levSetConfiguration);
    % fprintf('done!\n')

    % Initialize springs and x and y potential functions
    L      = zeros(g.shape);
    intVol = zeros(g.shape);

    MAX_ALLOWABLE_GRAD_PHI = 3;
    tNow  = 0;
    small = 100 * eps;

    %----------------------------------------------------------------------
    dtstp=1;
    dt = 0.001;
    % Start running the simulation
    tf = size(ext_Signal_norm,1)*dt;
    for stepi = 1:1:round((tf/dt))
        % fprintf('\nStep %g ',stepi)

        tStart = cputime;
        tSpan = [tNow,tNow+dt];

        % Data to be dealt with during the simulation
        ENData = [membranePotentialFunction(:) L(:)];

        % ExternalSignal
        schemeData.externalSignal = ext_Signal_norm(floor(stepi/dtstp),:);

        % Run the simulation
        [tNext, ENData] = feval(integratorFunc, @termChemoattractant,...
            tSpan, ENData, integratorOptions, schemeData);

        % Reinit membraneData if necessary
        phi = reshape(ENData(:,1), g.shape);
        L   = reshape(ENData(:,2), g.shape);
        myNormGradPhi = computeNormOfGradPhi(g, phi, schemeData, -1);
        myNorm = max(myNormGradPhi(:));
        if myNorm > 20
            % When spikes happen, re-calculate phi to be the signed distance
            % function based only on the membrane.
            % fprintf('\nNorm of grad phi is %g, recreating signed distance function...', myNorm);
            membranePoints = contourc(g.xs{1}(:,1),g.xs{2}(1,:),phi',[0 0]);
            membranePoints = membranePoints(:,2:(membranePoints(2,1)));
            phi = createSignedDistanceFunction(levSetConfiguration, ...
                membranePoints);
            % fprintf('Done!');
            myNormGradPhi = computeNormOfGradPhi(g, phi, schemeData, -1);
            myNorm = max(myNormGradPhi(:));
        end
        while(myNorm > MAX_ALLOWABLE_GRAD_PHI)
            % fprintf('\nNorm of grad phi is %g, calling Reinit...', myNorm);
            [phi, tElapsed] = reinitOnTheFly(g, phi, 0.2, 'low');
            myNormGradPhi = computeNormOfGradPhi(g, phi, schemeData, -1);
            myNorm = max(myNormGradPhi(:));
            % fprintf('Done!');
        end

        % Get back the correctly shaped data array
        membranePotentialFunction = phi;
        % Update membranePoints with new boundary
        membranePoints = contourc(g.xs{1}(:,1),g.xs{2}(1,:), ...
            membranePotentialFunction', [0 0]);

        % Pulls off the first 0 contour (in case there are multiple ones...
        %   ie. a cell splitting...)
        membranePoints = membranePoints(:,2:(membranePoints(2,1)));
        if (polygon_area_2d(membranePoints)>=0)
            %if curve is ccw, make it clockwise
            membranePoints = membranePoints(:,end:-1:1);
        end
        % Update the spring potential function
        memSpring = interp2(g.xs{1}',g.xs{2}', L',membranePoints(1,:), membranePoints(2,:));

        %  Filter out the spring state larger than 5 or smaller than -5
        memSpring=min(5,max(memSpring,-5));
        spring = griddata(membranePoints(1,:), membranePoints(2,:), ...
            memSpring, g.xs{1}, g.xs{2}, 'nearest');
        % Update the time step
        tNow = tNext;
        % Variables for length of protrusion versus time
        [geom, iner, cpmo] = polygeom(membranePoints(1,:)',membranePoints(2,:)');

        if ((dtstp==1)||(mod(stepi,dtstp)==1))
            savefile = fullfile(folder_name,strcat(file_type,'_lsm_output.mat'));
            save(savefile,'schemeData');
        end
        %
        n = schemeData.n;
        schemeData.membranePoints{n} = membranePoints;
        schemeData.memSpring{n} = memSpring;
        schemeData.surfArea(n)  = geom(1);
        schemeData.perimeter(n) = geom(4);
        schemeData.centerX(n)   = geom(2);
        schemeData.centerY(n)   = geom(3);
        schemeData.time(n)      = tNow;
        schemeData.n            = n+1;
        % Finish this time step
        fprintf('\n Step %g/%g, computation time %g seconds',...
            stepi, round((tf/dt)),cputime-tStart);
    end

end

    
    