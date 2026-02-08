%% function run_for_OE_kymo
% generates kymographs for Overexpression phenotype


    clearvars 
    close all
    clc
    warning('off','all');


    rng(200)

%% Spatial domain
Rc   = 5;   % radius of cell (um)
p.Np = 314; % number of spatial point along the perimeter
p.x = linspace(0,2*pi*Rc,p.Np); 
p.dx = mean(diff(p.x));   % spatial steps (um)



p.scale_spatial = 0.2424;

%% Time
p.scaling = 0.95;%1*0.9;
p.scaling2 = 1;
p.t0 = 0;   % Length of simulation (s)
p.tf = 900;   % Length of simulation (s)
p.dt = 0.001; % s
t = p.t0:p.dt:p.tf;

 cell_type_list = {'OE'};
 tic
    for k = 1
         p.cell_type = cell_type_list{k};

%% Parameters
%--------------------------------------------------------------------------
% dRas_dt
%--------------------------------------------------------------------------
    p.a1 = 0.133*p.scaling*p.scaling2;
    p.a2 = 2.85*p.scaling*p.scaling2;
    p.a3 = 1.14*p.scaling*p.scaling2;
    p.a4 = 1.05*440; 
    p.a5 = 0.0010*p.scaling*p.scaling2;
%--------------------------------------------------------------------------
% dPIP2_dt
%--------------------------------------------------------------------------  
    p.b1 = 0.0144*p.scaling*p.scaling2;
    p.b2 = 2134*p.scaling*p.scaling2;
    p.b3 = 1.7493*p.scaling*p.scaling2; 
    p.w  = 1.3*4.2*p.scaling*p.scaling2; % activation of PIP2 by PIP5K (membrane)
%--------------------------------------------------------------------------
% dPKB_dt 
%--------------------------------------------------------------------------
    p.c1 = 0.0744*p.scaling*p.scaling2;
    p.c2 = 0.8640*p.scaling*p.scaling2;
%--------------------------------------------------------------------------
% dPIP5K_dt 
%--------------------------------------------------------------------------
    p.d1 = 1*p.scaling*p.scaling2;
    p.d2 = 0.02;
    p.d3 = 7.5*0.05*p.scaling*p.scaling2;
%--------------------------------------------------------------------------
% dActin_dt 
%--------------------------------------------------------------------------
    p.ac1 = 1*0.05*p.scaling*p.scaling2;
    p.ac2 = 1*0.05*p.scaling*p.scaling2;
%--------------------------------------------------------------------------
% dMyosin_dt 
%--------------------------------------------------------------------------
    p.m1  = 0.8*0.05*p.scaling*p.scaling2; 
    p.m2  = 0.8*0.05*p.scaling*p.scaling2;
%--------------------------------------------------------------------------
% Diffusion coefficients
%--------------------------------------------------------------------------
    p.DRas  = 1*1.1*0.03*p.scaling*p.scale_spatial*p.scaling2;
    p.DPIP2 = 1*0.03*p.scaling*p.scale_spatial*p.scaling2;
    p.DPKB  = 3*0.03*p.scaling*p.scale_spatial*p.scaling2;

    p.DActin  = 2*0.0025*p.scaling*p.scale_spatial*p.scaling2;
    p.DMyosin = 2*2*0.0025*p.scaling*p.scale_spatial*p.scaling2;

    p.DPIP5K_cyto = 20*0.0025*400*p.scaling*p.scale_spatial*p.scaling2;
    p.DPIP5K_mem  = 10*10*0.0025*p.scaling*p.scale_spatial*p.scaling2;
%--------------------------------------------------------------------------
% Noise p.scaling factor
%--------------------------------------------------------------------------
    p.alpha = 1/12*(1/0.9)*p.scaling;

   
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    switch p.cell_type
       
         case 'OE'
            p.a = 0.05*0.25*10/10; % Ras activation by Actin: p.a; change if needed
            p.m = 0.05*0.25*6/10;  % Ras inhibition by Myosin: p.m; change if needed
            p.c2 = 1.1*0.8640*p.scaling*p.scaling2;
 
    end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% nullclines and equilibrium values of Ras, PIP2, PKB, Actin, Myosin
%--------------------------------------------------------------------------
    p = nullclines(p); 
%--------------------------------------------------------------------------
% initial values for PIP5K on membrane and in cytosol
%--------------------------------------------------------------------------
    switch p.cell_type
            
            case 'OE'
                % higher initial amount of PIP5K
                p.PIP5K_total = 1.7*10/126;
                p.PIP5K_mem0  = 1.7*9/126;
                
            
    end




    p.PIP5K_cyto0 = p.PIP5K_total - p.PIP5K_mem0;
%--------------------------------------------------------------------------

    p.Tmem0 = 0;

    initvalues = [p.Ras0, p.PIP20, p.PKB0, p.PIP5K_mem0, p.PIP5K_cyto0, p.Actin0, p.Myosin0, p.Tmem0]; 
    

%% SDE toolbox

    
    timeRange = t;
    problem = 'kymo'; % name of the sde file
    numsim  = p.Np; % number of spatial points
    sdetype = 'Ito';
    numdepvars = length(initvalues);

   

       
    Output = SDE_euler_deb(initvalues,problem,timeRange,numdepvars,numsim,sdetype,p);
    

% Storing the data 

        T_sampling = 1/p.dt/10;


        Ras    = Output(1:T_sampling:end,1:numdepvars:end);
        PIP2   = Output(1:T_sampling:end,2:numdepvars:end);
        PKB    = Output(1:T_sampling:end,3:numdepvars:end);


        PIP5K_mem    = Output(1:T_sampling:end,4:numdepvars:end);
        PIP5K_cyto   = Output(1:T_sampling:end,5:numdepvars:end);

        Actin    = Output(1:T_sampling:end,6:numdepvars:end);
        Myosin   = Output(1:T_sampling:end,7:numdepvars:end);

        Tmem   = Output(1:T_sampling:end,8:numdepvars:end);



        T = t(1:T_sampling:end);
%% Saving the data
        % target_folder = '../Kymograph_Data';
        % if ~exist(target_folder,"dir")
        %     mkdir(target_folder)
        % end
        % filename = fullfile(target_folder,strcat(p.cell_type,'_data.mat'));
        % 
        % save(filename,"p","Ras","PIP2","PKB","PIP5K_mem","PIP5K_cyto","Actin","Myosin","T");

    end
        Elapsed_time =  toc;

        fprintf('Elapsed time in min = %f\n',Elapsed_time/60)

%% Plot

d = 100;
close all
T_initial = 0;

figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(smoothdata(Ras)',d,1))

colorbar
clim([0 0.2])
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 (p.tf - T_initial)])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('Time','fontsize',20)
ylabel('Cell perimeter','fontsize',20)
title([cell_type_list{k},', Ras Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end

%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(smoothdata(PIP2,"gaussian")',d,1))
colorbar
clim([0 1])
mycolormap = magma(100);
colormap(mycolormap)
% colormap(flipud(sky))
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', PIP2 Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(PIP2)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(PKB',d,1))
colorbar
clim([0 1])
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', PKB Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(PKB)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(Actin',d,1))
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', Actin Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
clim([0 0.7])
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(Myosin',d,1))
% surf(T-T_initial,p.x,circshift(Myosin',d,1),'edgecolor','interp')
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', Myosin Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(Tmem',d,1))
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', Tmem Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(Actin'-Tmem',d,1))
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', Actin - Tmem Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(2*Actin'-Myosin',d,1))
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
% title([cell_type_list{k},', Actin - Myosin Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
clim([-0.25 1.25])
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
sum_Ras = sum(Ras,2);
sum_PIP5K_mem = sum(PIP5K_mem,2);
sum_PIP5K_cyto = sum(PIP5K_cyto,2);
total_PIP5K = sum_PIP5K_mem + sum_PIP5K_cyto;
% plot(T-T_initial,smooth(sum_Ras),'-','LineWidth',2,'Color',[36 136 36]/255)
h(1) = plot(T-T_initial,total_PIP5K,'-','LineWidth',3,'Color',0.7*[1 1 1]);
h(2) = plot(T-T_initial,sum_PIP5K_mem,'-','LineWidth',2,'Color',[220, 20, 60]/255);
h(3) = plot(T-T_initial,sum_PIP5K_cyto,'-.','LineWidth',2,'Color',[0, 150, 255]/255);
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
title([cell_type_list{k},', Total Amount of PIP5K'],'fontsize',20)
legend(h,{'Total PIP5K','PIP5K mem','PIP5K cyto'},'fontsize',16,'EdgeColor','none')
xlim([0 p.tf - T_initial])
ylim([0 40])

%%
ratio_PIP5K_mem2cyto = (mean(sum_PIP5K_mem)/mean(sum_PIP5K_cyto));
fprintf('PIP5K ratio membrane to cytosol = %0.2f\n',round(ratio_PIP5K_mem2cyto,2))

%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(PIP5K_mem',d,1))
% surf(T-T_initial,p.x,circshift(Myosin',d,1),'edgecolor','interp')
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', PIP5K (membrane) Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
% clim([0 max(PIP5K_mem(:) + PIP5K_cyto(:))])
clim([0 0.12])
%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(PIP5K_cyto',d,1))
% surf(T-T_initial,p.x,circshift(Myosin',d,1),'edgecolor','interp')
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', PIP5K (cytosol) Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
% clim([0 max(PIP5K_mem(:) + PIP5K_cyto(:))])


%%
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',16)
hold on
imagesc(T-T_initial,p.x,circshift(PIP5K_mem' + PIP5K_cyto',d,1))
% surf(T-T_initial,p.x,circshift(Myosin',d,1),'edgecolor','interp')
colorbar
colormap(gray)
mycolormap = magma(100);
colormap(mycolormap)
xlim([0 p.tf - T_initial])
ylim([0 p.x(end)])
xticks(0:300:p.tf);
xticklabels(num2cell(0:300:p.tf))
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title([cell_type_list{k},', PIP5K Total Kymograph'],'fontsize',20,'fontweight','n')
for kk = 1:60:length(Ras)
   plot([kk+.1,kk+.1],[0,p.x(end)],'y-');
end
clim([0 max(PIP5K_mem(:) + PIP5K_cyto(:))])