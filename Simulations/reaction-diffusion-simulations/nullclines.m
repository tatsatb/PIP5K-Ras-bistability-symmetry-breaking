function newp = nullclines(p)

syms Ras 
newp = p;

%--------------------------------------------------------------------------
% quasi-steady state values
%--------------------------------------------------------------------------
PIP2  = p.b3/(p.b1 + p.b2*Ras);
PKB  = p.c2*Ras/p.c1;
Actin  = p.ac2*PKB/p.ac1;
Myosin = p.m2*PIP2/p.m1;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Reaction fluxes for Ras
%--------------------------------------------------------------------------
J1 = (p.a1 + p.a2*PKB)*Ras;
J2 = p.a3./((p.a4^2*PIP2.^2 + 1)) + p.a5; 
J1 = J1*(1 + p.m*Myosin); % modification by Myosin inhibition
J2 = J2*(1 + p.a*Actin);  % modification by Actin activation


eqn = -J1+J2==0;
S0 = vpasolve(eqn,Ras);
S0 = double(S0);
S0 = max(S0(S0==real(S0)));

newp.Ras0 = S0(S0>0);
newp.PIP20 = p.b3/(p.b1 + p.b2*newp.Ras0);
newp.PKB0 = p.c2*newp.Ras0/p.c1;
newp.Actin0 = p.ac2*newp.PKB0/p.ac1;
newp.Myosin0 = p.m2*newp.PIP20/p.m1;

% plot the nullclines
    Ras  = [1e-5:1e-5:0.02 0.021:1e-3:5];
    PKB = (p.c2*Ras)/p.c1;
    PIP2 = p.b3./(p.b1 + p.b2*Ras);
    Actin = p.ac2*PKB/p.ac1;
    Myosin = p.m2*PIP2/p.m1;

    J2 = p.a3./((p.a4^2*PIP2.^2 + 1)) + p.a5;
    J2 = J2.*(1 + p.a*Actin);
    f = (J2./(1 + p.m*Myosin)./Ras - p.a1)/p.a2;
    r = p.c2*Ras/p.c1;

    figure('Color','white')
    hold on
    plot(Ras,f,'Color',[36 136 36]/255,'LineWidth',2)
    plot(Ras,r,'Color',[0 0 255]/255,'LineWidth',2)
    plot(newp.Ras0,newp.PKB0,'o','MarkerSize',10,'MarkerFaceColor',[0 1 1],'color','k','LineWidth',2)
    set(gca,'Yscale','log')
    set(gca,'Xscale','log')
    xlabel('Ras','FontSize',16)
    ylabel('PKB','FontSize',16)
end % funnction p=nullclines(p)