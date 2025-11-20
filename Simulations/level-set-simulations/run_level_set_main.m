%%function run_level_set_main
% run this code for level set simulation of the kymograph data stored in
% "Kymograph_Data"


    clearvars
    close all
    clc
    warning('off','all')


    param = param_list;




%% change the file type

    file_type = 'WT';
    % file_type = 'pikl_fan';
    % file_type = 'pikl_oscillator';
    % file_type = 'OE';

    folder_path_to_kymodata = '../Kymograph_Data';



    file_list = dir(folder_path_to_kymodata);

    for i = 1:length(file_list)

        if contains(file_list(i).name,file_type)
            db_data = load(fullfile(folder_path_to_kymodata,file_list(i).name));
        end

    end
%--------------------------------------------------------------------------
% extract the segment of data for Level Set Simulation
%--------------------------------------------------------------------------
    switch file_type
        case 'WT'
            idx1 = 1200;
            idx2 = 4200;

         case 'pikl_fan'           
            idx1 = 4200;
            idx2 = 7200;

% fan cells covers longer spatial range so update accordingly
            param.xleft_grid = -50;
            param.xright_grid = 10;
            param.ybottom_grid = -20;
            param.ytop_grid = 20;
           
         case 'pikl_oscillator'            
            idx1 = 1050;
            idx2 = 4050;
            param.volCKp = 0.1*param.volCKp;            

         case 'OE'
            idx1 = 1200;
            idx2 = 4200;
       
     end
%--------------------------------------------------------------------------
% use PKB data as input to the Level Set
%--------------------------------------------------------------------------
    input = db_data.PKB;
    input = input(idx1:idx2,:);
    input = circshift(input,150,2);

% Scaling of the signal    
    ext_Signal = input;
    T_int = 10;
    ext_Signal  = ext_Signal(1:T_int:end,:);

    ext_Signal_norm = min(max((ext_Signal - 0.3),0),0.25);

% Plots to check the scaled input to Level Set
    plot_time2 = 0:1:size(ext_Signal_norm,1)-1;
    plot_time1 = (0:1:size(input,1)-1)*0.1;
    plot_x = 1:1:size(ext_Signal_norm,2);
    figure('Color','white')
    subplot(211)
    set(gca,'linewidth',1.5,'fontsize',12)
    hold on
    imagesc(plot_time1,plot_x,  input')
    colormap jet
    colorbar
    hold on;
    xlabel('Time (s)');
    ylabel('Perimeter Points');
    title_str = 'PKB';
    title(title_str,'FontWeight','normal');
    axis tight
     
    subplot(212)
    set(gca,'linewidth',1.5,'fontsize',12)
    hold on
    imagesc(plot_time2,plot_x,ext_Signal_norm')
    colormap jet
    colorbar
    hold on;
    % clim([-0.5 1])
    xlabel('Time (s)');
    ylabel('Perimeter Points');
    title('Input to LS','FontWeight','normal');
    axis tight

%%

% cd ./LSMCode
addpath("LSM_Toolbox/")
tic
folder_name = 'LSM_output';
LSM_moving_cell(param,file_type,ext_Signal_norm,folder_name)
elapsed_time = toc;
fprintf('Elapsed time in hr = %f\n',elapsed_time/3600);


% end % end of function