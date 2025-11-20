function param = param_list

%% ---------------------------------------------------------------------
% LSM parameters:
% The cell's mechanical properties; from Yang et al., BMC Sys Bio, 2008
%---------------------------------------------------------------------
    Rc  = 5;   % Initial cell radius in um
    km   = 0.098; % s*nN/um^2
    gamm = 0.064; % s*nN/um^3
    gamc = 6.09;  % s*nN/um^3

    volCKp = 20.0; % (dimensionless);
    volCKi = 1.0;  % (dimensionless)
    gamTen = 2; % nN/um; from Yang et al., BMC Sys Bio, 2008
    gs = 11; % Grid numbers per um 
%---------------------------------------------------------------------
%% Grid size for LSM
%---------------------------------------------------------------------
    xleft_grid   = -35;
    xright_grid  = +35;
    ybottom_grid = -35;
    ytop_grid    = +35;

%% ---------------------------------------------------------------------


param.Rc = Rc;
param.km = km;
param.gamm = gamm;
param.gamc = gamc;
param.volCKp = volCKp;
param.volCKi = volCKi;
param.gamTen = gamTen;
param.gs = gs;

param.xleft_grid = xleft_grid;
param.xright_grid = xright_grid;
param.ybottom_grid = ybottom_grid;
param.ytop_grid = ytop_grid;
end

