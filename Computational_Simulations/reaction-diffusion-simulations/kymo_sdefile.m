

function [out1,out2,out3] = kymo_sdefile(t,x,flag,initvalues,SDETYPE,NUMDEPVARS,NUMSIM,param)


% Initial condition


    Ras_zero    = initvalues(1);     
    PIP2_zero   = initvalues(2);    
    PKB_zero    = initvalues(3);     

    PIP5K_mem_zero   = initvalues(4);     
    PIP5K_cyto_zero   = initvalues(5);     

    Actin_zero      = initvalues(6);  
    Myosin_zero     = initvalues(7); 
    Tmem_zero       = initvalues(8);  


  


if nargin < 3 || isempty(flag)
    

    xsplitted = cell(1,NUMDEPVARS);
    for i=1:NUMDEPVARS
        xsplitted{i} = x(i:NUMDEPVARS:end);
    end
    

  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  %:::::::::::::::::::::::::::::::  DEFINE HERE THE SDE  :::::::::::::::::::::::::::::
  %::::::::::::: (define the initial conditions at the bottom of the page) :::::::::::
  %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  Ras      = xsplitted{1}; 
  PIP2     = xsplitted{2};
  PKB      = xsplitted{3};
  
  PIP5K_mem  = xsplitted{4}; 
  PIP5K_cyto = xsplitted{5}; 

  Actin    = xsplitted{6};
  Myosin   = xsplitted{7}; 
  Tmem     = xsplitted{8}; 

 
           
   switch upper(SDETYPE)
   case 'ITO'
      

   
   


%=======================================================       

  
  
 
  
%--------------------------------------------------------------------------
    J1_Ras = (param.a1 + param.a2*PKB).*(1 + param.m*Myosin).*Ras;
    % J2_Ras = (param.a3./((param.a4.^2*PIP2.^2 + 1)) + param.a5).*(1 + param.a*(Actin-param.Tmem*Tmem));
    J2_Ras = (param.a3./((param.a4.^2*PIP2.^2 + 1)) + param.a5).*(1 + param.a*Actin - param.Tmem*Tmem);
    
    drift_Ras  = - J1_Ras + J2_Ras + param.DRas./(param.dx^2)*nabla_SDE(Ras);
    noise_Ras  = sqrt(J1_Ras + J2_Ras)*param.alpha;
%--------------------------------------------------------------------------

    J1_PIP2 = (param.b1 + param.b2*Ras).*PIP2;
    J2_PIP2 = param.b3.*(1 + param.b4*PIP5K_mem);    
    
    drift_PIP2  = - J1_PIP2 + J2_PIP2 + param.DPIP2./(param.dx^2)*nabla_SDE(PIP2);
    noise_PIP2  = sqrt(J1_PIP2 + J2_PIP2)*param.alpha;
%--------------------------------------------------------------------------
    J1_PKB = param.c1*PKB; 
    J2_PKB = param.c2*Ras;

    drift_PKB = - J1_PKB + J2_PKB + param.DPKB./(param.dx^2)*nabla_SDE(PKB);
    noise_PKB  = sqrt(J1_PKB + J2_PKB)*param.alpha;
%--------------------------------------------------------------------------

    J1_PIP5K = param.d1*PIP5K_mem.*Ras.^2./(param.d2^2+Ras.^2); %Ras(F) releases PI5K(mem)
    J2_PIP5K = param.d3*PIP5K_cyto;


    drift_PIP5K_mem  = - J1_PIP5K + J2_PIP5K + param.DPIP5K_mem./(param.dx^2).*nabla_SDE(PIP5K_mem);
    drift_PIP5K_cyto =   J1_PIP5K - J2_PIP5K + param.DPIP5K_cyto./(param.dx^2).*nabla_SDE(PIP5K_cyto);
    
    noise_PI5K_mem   = 0;
    noise_PI5K_cyto  = 0;

%--------------------------------------------------------------------------

    J1_Actin  = param.ac1*Actin;
    J2_Actin  = param.ac2*PKB;
    drift_Actin = - J1_Actin + J2_Actin + param.DActin./(param.dx^2)*nabla_SDE(Actin);
    noise_Actin  = 1*sqrt(J1_Actin + J2_Actin)*param.alpha;
%--------------------------------------------------------------------------
    J1_Myosin  = param.m1*Myosin;
    J2_Myosin  = param.m2*PIP2;
    drift_Myosin = - J1_Myosin + J2_Myosin + param.DMyosin./(param.dx^2)*nabla_SDE(Myosin);
    noise_Myosin  = 1*sqrt(J1_Myosin + J2_Myosin)*param.alpha;
%--------------------------------------------------------------------------
    J1_Tmem  = 0.125*Tmem;
    J2_Tmem  = 0.125*mean(PKB);
    drift_Tmem = - J1_Tmem + J2_Tmem;% p.DTmem./(p.dx^2)*nabla_SDE(Tmem);
    noise_Tmem  = 0*sqrt(J1_Tmem + J2_Tmem)*param.alpha;
%--------------------------------------------------------------------------   
   end
   
    out1 = zeros(1,NUMDEPVARS*NUMSIM);
    out1(1:NUMDEPVARS:end) = drift_Ras;
    out1(2:NUMDEPVARS:end) = drift_PIP2;
    out1(3:NUMDEPVARS:end) = drift_PKB;
  
    out1(4:NUMDEPVARS:end) = drift_PIP5K_mem;
    out1(5:NUMDEPVARS:end) = drift_PIP5K_cyto;
    
    out1(6:NUMDEPVARS:end) = drift_Actin;   
    out1(7:NUMDEPVARS:end) = drift_Myosin;

    out1(8:NUMDEPVARS:end) = drift_Tmem;

%--------------------------------------------------------------------------    
    out2 = zeros(1,NUMDEPVARS*NUMSIM);
    out2(1:NUMDEPVARS:end) = noise_Ras;
    out2(2:NUMDEPVARS:end) = noise_PIP2;
    out2(3:NUMDEPVARS:end) = noise_PKB;
   
    out2(4:NUMDEPVARS:end) = noise_PI5K_mem;
    out2(5:NUMDEPVARS:end) = noise_PI5K_cyto;

    out2(6:NUMDEPVARS:end) = noise_Actin;
    out2(7:NUMDEPVARS:end) = noise_Myosin;

    out2(8:NUMDEPVARS:end) = noise_Tmem;
%--------------------------------------------------------------------------        
    out3 = zeros(1,NUMDEPVARS*NUMSIM);
   
   %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
else
    
    switch(flag)
    case 'init'  
        out1 = t;
        
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%::::::::::::::::::::::  DEFINE HERE THE SDE INITAL CONDITIONS  :::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

out2 = [Ras_zero PIP2_zero PKB_zero PIP5K_mem_zero PIP5K_cyto_zero Actin_zero Myosin_zero Tmem_zero];   % write here the SDE initial condition(s)
        
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: :::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        out3 = [];
        
        
    otherwise
        error(['Unknown flag ''' flag '''.']);
    end
end
