%% function run_for_WT
% generates kymographs for WT phenotype

    clearvars 
    close all
    clc
    warning('off','all');

    rng(100)

%% load parameters

param = param_list;

%% Simulation Time

t = param.t0:param.dt:param.tf;

 cell_type_list = {'WT'};
 tic
    for k = 1
         param.cell_type = cell_type_list{k};

%--------------------------------------------------------------------------

    switch param.cell_type
        case 'WT'
            param.a = 0.2750; % Ras activation by Actin: p.a; change if needed
            param.m = 0.1500;  % Ras inhibition by Myosin: p.m; change if needed
 
    end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% nullclines and equilibrium values of Ras, PIP2, PKB, Actin, Myosin
%--------------------------------------------------------------------------
    param = nullclines(param); 
%--------------------------------------------------------------------------
% initial values for PIP5K on membrane and in cytosol
%--------------------------------------------------------------------------
    switch param.cell_type
            case 'WT'
                param.PIP5K_total = 0.0794;
                param.PIP5K_mem0 = 0.0714; % 90% PIP5K are on membrane initially
           
    end




    param.PIP5K_cyto0 = param.PIP5K_total - param.PIP5K_mem0;
%--------------------------------------------------------------------------

    param.Tmem0 = 0;

    initvalues = [param.Ras0, param.PIP20, param.PKB0, param.PIP5K_mem0, param.PIP5K_cyto0, param.Actin0, param.Myosin0, param.Tmem0]; 
    

%% SDE toolbox

    
    timeRange = t;
    problem = 'kymo'; % name of the sde file
    numsim  = param.Np; % number of spatial points
    sdetype = 'Ito';
    numdepvars = length(initvalues);
   
    Output = SDE_euler_deb(initvalues,problem,timeRange,numdepvars,numsim,sdetype,param);
    
%% Store sampled data

    T_sampling = round(1/(10*param.dt));   % sample every 0.1 time units

    speciesNames = {'Ras','PIP2','PKB','PIP5K_mem','PIP5K_cyto','Actin','Myosin','Tmem'};

    for i = 1:numel(speciesNames)

        S.(speciesNames{i}) = Output(1:T_sampling:end, i:numdepvars:end);

    end

    Ras        = S.Ras;
    PIP2       = S.PIP2;
    PKB        = S.PKB;
    PIP5K_mem  = S.PIP5K_mem;
    PIP5K_cyto = S.PIP5K_cyto;
    Actin      = S.Actin;
    Myosin     = S.Myosin;
    Tmem       = S.Tmem;

    T = t(1:T_sampling:end);

%% Save kymograph data

        target_folder = '../Kymograph_Data_new';

        if ~exist(target_folder,'dir')
            mkdir(target_folder);
        end

        filename = fullfile(target_folder, sprintf('%s_data.mat', param.cell_type));

        save(filename, ...
            'param','T', ...
            'Ras','PIP2','PKB', ...
            'PIP5K_mem','PIP5K_cyto', ...
            'Actin','Myosin','Tmem');
    end
        Elapsed_time =  toc;

        fprintf('Elapsed time in min = %f\n',Elapsed_time/60)

%% Plot


d = 100;
T_initial = 0;
close all

time = T - T_initial;
xpos = param.x;
tf_plot = param.tf - T_initial;
xtick_vals = 0:300:param.tf;
event_step = 60;
event_color = 'y-';

commonFont = 16;
labelFont = 20;
lw_axis = 1;

cmap = magma(100);

%% Helper function for kymographs
plot_kymo = @(data,title_str,clim_vals) local_kymograph( ...
    time, xpos, data, d, tf_plot, xtick_vals, ...
    title_str, clim_vals, cmap, event_step, event_color, ...
    commonFont, labelFont, lw_axis);

%% Kymographs
plot_kymo(smoothdata(Ras)', ...
    [cell_type_list{k}, ', Ras Kymograph'], [0 0.2]);

plot_kymo(smoothdata(PIP2,"gaussian")', ...
    [cell_type_list{k}, ', PIP2 Kymograph'], [0 1]);

plot_kymo(PKB', ...
    [cell_type_list{k}, ', PKB Kymograph'], [0 1]);

plot_kymo(Actin', ...
    [cell_type_list{k}, ', Actin Kymograph'], [0 0.7]);

plot_kymo(Myosin', ...
    [cell_type_list{k}, ', Myosin Kymograph'], []);

plot_kymo(Tmem', ...
    [cell_type_list{k}, ', Tmem Kymograph'], []);

plot_kymo(Actin' - Tmem', ...
    [cell_type_list{k}, ', Actin - Tmem Kymograph'], []);

plot_kymo(2*Actin' - Myosin', ...
    [cell_type_list{k}, ', 2 Actin - Myosin Kymograph'], [-0.25 1.25]);

%% Total PIP5K time series
sum_Ras        = sum(Ras,2);
sum_PIP5K_mem  = sum(PIP5K_mem,2);
sum_PIP5K_cyto = sum(PIP5K_cyto,2);
total_PIP5K    = sum_PIP5K_mem + sum_PIP5K_cyto;

figure('Color','white')
set(gca,'FontWeight','normal','LineWidth',lw_axis,'FontSize',commonFont)
hold on

h(1) = plot(time,total_PIP5K,'-','LineWidth',3,'Color',0.7*[1 1 1]);
h(2) = plot(time,sum_PIP5K_mem,'-','LineWidth',2,'Color',[220 20 60]/255);
h(3) = plot(time,sum_PIP5K_cyto,'-.','LineWidth',2,'Color',[0 150 255]/255);

xlim([0 tf_plot])
ylim([0 40])
xticks(xtick_vals)
xticklabels(num2cell(xtick_vals))

xlabel('Time','FontSize',labelFont)
title([cell_type_list{k}, ', Total Amount of PIP5K'], ...
    'FontSize',labelFont,'FontWeight','normal')

legend(h,{'Total PIP5K','PIP5K mem','PIP5K cyto'}, ...
    'FontSize',commonFont,'EdgeColor','none')

%% PIP5K membrane-to-cytosol ratio
ratio_PIP5K_mem2cyto = mean(sum_PIP5K_mem) / mean(sum_PIP5K_cyto);
fprintf('PIP5K ratio membrane to cytosol = %0.2f\n', ratio_PIP5K_mem2cyto)

%% PIP5K kymographs
plot_kymo(PIP5K_mem', ...
    [cell_type_list{k}, ', PIP5K (membrane) Kymograph'], [0 0.12]);

plot_kymo(PIP5K_cyto', ...
    [cell_type_list{k}, ', PIP5K (cytosol) Kymograph'], []);

plot_kymo(PIP5K_mem' + PIP5K_cyto', ...
    [cell_type_list{k}, ', PIP5K Total Kymograph'], ...
    [0 max(PIP5K_mem(:) + PIP5K_cyto(:))]);


%% Local function
function local_kymograph(time,xpos,data,d,tf_plot,xtick_vals, ...
    title_str,clim_vals,cmap,event_step,event_color, ...
    commonFont,labelFont,lw_axis)

    figure('Color','white')
    set(gca,'FontWeight','normal','LineWidth',lw_axis,'FontSize',commonFont)
    hold on

    imagesc(time,xpos,circshift(data,d,1))
    axis xy

    colorbar
    colormap(cmap)

    if ~isempty(clim_vals)
        clim(clim_vals)
    end

    xlim([0 tf_plot])
    ylim([0 xpos(end)])
    xticks(xtick_vals)
    xticklabels(num2cell(xtick_vals))

    xlabel('Time','FontSize',labelFont)
    ylabel('Cell perimeter','FontSize',labelFont)
    title(title_str,'FontSize',labelFont,'FontWeight','normal')

    for kk = 1:event_step:length(time)
        plot([kk kk],[0 xpos(end)],event_color,'LineWidth',1)
    end
end
