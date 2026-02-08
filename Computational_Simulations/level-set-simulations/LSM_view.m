%% function LSM_view
% Script to visualize level set simulations

clearvars
close all
clc



addpath('../Kymograph_Code/')

%% Uncomment accordingly

file_type = 'WT';
% file_type = 'pikl_fan';
% file_type = 'pikl_oscillator';
% file_type = 'OE';

%--------------------------------------------------------------------------
% load level set output and kymograph data and name the video output file:
%--------------------------------------------------------------------------
% folder_path_to_lsmdata  = '../Level_Set_Data';
% 
%     file_list_lsm = dir(folder_path_to_lsmdata);
% 
%     for i = 1:length(file_list_lsm)
% 
%         if contains(file_list_lsm(i).name,file_type)
%             load(fullfile(folder_path_to_lsmdata,file_list_lsm(i).name));
%         end
% 
%     end
load('LSM_output/WT_lsm_output.mat') 

%--------------------------------------------------------------------------
folder_path_to_kymodata = '../Kymograph_Data';        
    file_list_kymo = dir(folder_path_to_kymodata);

    for i = 1:length(file_list_kymo)

        if contains(file_list_kymo(i).name,file_type)
            db_data = load(fullfile(folder_path_to_kymodata,file_list_kymo(i).name));
        end

    end
%--------------------------------------------------------------------------
    
    video_folder_name = 'LSM_video';

    if ~exist(video_folder_name,'dir')
        
        mkdir(video_folder_name)

    end
    video_file_string = fullfile(video_folder_name,strcat(file_type,'_video'));
%--------------------------------------------------------------------------
    switch file_type
        case 'WT'
            idx1 = 1200;
            idx2 = 4200;

         case 'pikl_fan'           
            idx1 = 4200;
            idx2 = 7200;

            param.xleft_grid = -50;
            param.xright_grid = 10;
            param.ybottom_grid = -20;
            param.ytop_grid = 20;
           
         case 'pikl_oscillator'            
            idx1 = 50;
            idx2 = 6050;       

         case 'OE'
            idx1 = 1200;
            idx2 = 4200;
       
     end
%--------------------------------------------------------------------------
    input = db_data.PKB;

    input = input(idx1:idx2,:);
    input = circshift(input,150,2);
    ext_Signal = input;
    delT = 10;
    ext_Signal  = ext_Signal(1:delT:end,:);
    Rc   = 5;   % radius of cell (um)
    p.Np = 314; % number of spatial point along the perimeter
    p.x = linspace(0,2*pi*Rc,p.Np);
    time = 0:1:(size(ext_Signal,1)-1);
%--------------------------------------------------------------------------
    
   
%%
% cmap = jet(length(schemeData.membranePoints)+50);


x_comb = cell(length(schemeData.membranePoints),1);
y_comb = cell(length(schemeData.membranePoints),1);
xc_comb = nan(length(schemeData.membranePoints),1);
yc_comb = nan(length(schemeData.membranePoints),1);

for i = 1:1:length(schemeData.membranePoints)
        x_comb{i} = schemeData.membranePoints{i}(1,:);
        y_comb{i} = schemeData.membranePoints{i}(2,:);
       
        polyin = polyshape(x_comb{i},y_comb{i});
        [xc_comb(i),yc_comb(i)] = centroid(polyin);

end


%%

T_init = 1;
cmap = jet(301+50);

vidWriter = VideoWriter(video_file_string, 'Motion JPEG AVI');
vidWriter.FrameRate = 10;
vidWriter.Quality = 100;
open(vidWriter);

for i = T_init:5:length(schemeData.membranePoints)
i  
        figure(1)
        set(gcf,'Visible', 'off');
        set(gcf,'color','black')
        % set(gcf, 'Position', [100, 100, 1000, 1000])
        subplot(2,2,1)
        hold on

        x = schemeData.membranePoints{i}(1,:);
        y = schemeData.membranePoints{i}(2,:);
        plot(x,y,'-','LineWidth',1.5,'Color','w');

        axis([-30 10 -10 20]);
        pbaspect([40 30 1])
        
        ax1 = gca;
        set(ax1, 'Color', 'k');            % Set axes background to black
        set(ax1, 'XColor', 'w');           % Set x-axis color to white
        set(ax1, 'YColor', 'w');           % Set y-axis color to white
        set(ax1, 'ZColor', 'w');           % For 3D plots
        xlabel('x position (\mum)','Color', 'w')
        ylabel('y position (\mum)','Color', 'w')
        sgtitle(['Time = ',num2str(i-1), ' s'],'FontWeight','normal','Color','w');
   
        subplot(2,2,2)
        hold on
        x = schemeData.membranePoints{i}(1,:);
        y = schemeData.membranePoints{i}(2,:);
        for j = 1:i
            plot(x_comb{j},y_comb{j},'-','LineWidth',1,'Color',[cmap(j+40, :), 0.5]);
            scatter(xc_comb(j), yc_comb(j), 15, cmap(j+40, :), 'filled','o','MarkerFaceAlpha',0.5);
        end
        plot(x,y,'-','LineWidth',2,'Color','w');
        
        axis([-30 10 -10 20]);
        pbaspect([40 30 1])
      
        ax2 = gca;
        set(ax2, 'Color', 'k');            % Set axes background to black
        set(ax2, 'XColor', 'w');           % Set x-axis color to white
        set(ax2, 'YColor', 'w');           % Set y-axis color to white
        set(ax2, 'ZColor', 'w');           % For 3D plots
        xlabel('x position (\mum)','Color', 'w')
        ylabel('y position (\mum)','Color', 'w')
       
        cb1 = colorbar;       
        colormap(ax2, cmap((T_init+40):341,:));
        hold on;
        % clim([T_init  length(schemeData.membranePoints)]-T_init )
        clim([T_init  301]-T_init )
            cb1.Color = 'w';                % Color of colorbar ticks and numbers
        cb1.Label.String = 'Time (s) ';     % Label text
        cb1.Label.Color = 'w';          % Label color
        cb1.Label.FontSize = 12;

        subplot(2,2,[3,4])
     
        imagesc((1:i)-1,p.x,flipud(ext_Signal(1:i,:)'))
        axis([0, 300, p.x(1) p.x(end)])
        
        cb2 = colorbar;
        hold on;
        clim([0 0.7])
        cb2.Color = 'w';                % Color of colorbar ticks and numbers
        cb2.Label.String = 'PKB Conc. (au) ';     % Label text
        cb2.Label.Color = 'w';          % Label color
        cb2.Label.FontSize = 12;

        xlabel('Time (s)');
        ylabel('Cell Perimeter (\mum)');
       
        ax3 = gca;
        set(ax3, 'Color', 'k');            % Set axes background to black
        set(ax3, 'XColor', 'w');           % Set x-axis color to white
        set(ax3, 'YColor', 'w');           % Set y-axis color to white
        set(ax3, 'ZColor', 'w');           % For 3D plots
        colormap(ax3, magma);
        % title(title_str,'FontWeight','normal');
        % return
    

  
     frame = getframe(gcf);
     writeVideo(vidWriter, frame);     % Write to output


  
    close all
end 
close all
close(vidWriter);
