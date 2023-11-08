% This function was written by the MIT 2.671 instructors

function [] = improvePlot()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot parameters
    % MATLAB treats mac and PC displays differently which can create
    % weird looking graphs. Here we handle system differences

    if ismac
        plot_width_in_px = 800;
        plot_height_in_px = 800;
        marker_size=15;
        marker_line_width=2.5;
        box_thickness = 3;
        axis_tick_font_size = 24;
        axis_label_font_size = 24;
        legend_font_size = 20;
        error_bar_cap_size = 15;
    else % (ispc || isunix)
        plot_width_in_px = 600;
        plot_height_in_px = 600;
        marker_size=10;
        marker_line_width=2.0;
        box_thickness = 2;
        axis_tick_font_size = 18;
        axis_label_font_size = 18;
        legend_font_size = 16;
        error_bar_cap_size = 10;
    end
    
    marker_outline = 'matching'; % could be 'black' or 'matching'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Use h as handle for current figure
    hFig = gcf;                    
    % Change figure background colour to white
    set(hFig, 'Color', 'white');

    % Make the figure bigger
    set(hFig, 'rend', 'painters', 'Units', 'pixels', 'pos', ...
        [100 100 plot_width_in_px plot_height_in_px]);

    % Grab the axes handle(s)
    axis_handles=findobj(hFig,'type','axe');

    % Iterate over all axes handle(s), this is useful if there are subplots
    for i = 1:length(axis_handles)
        ax = axis_handles(i);

        % Change default font size (tick labels, legend, etc.)
        set(ax, 'FontSize', axis_tick_font_size, 'FontName', 'Arial', 'LineWidth', box_thickness);
        
        set(ax, 'Box', 'on');

        % Change font size for axis text labels
        set(get(ax, 'XLabel'),'FontSize', axis_label_font_size, 'FontWeight', 'Bold');
        set(get(ax, 'YLabel'),'FontSize', axis_label_font_size, 'FontWeight', 'Bold');
        
        try % try statement to avoid error with categorical axes
        ax.XRuler.Exponent = 0; % Remove exponential notation from the X axis
        ax.YRuler.Exponent = 0; % Remove exponential notation from the Y axis
        catch
        end
        
    end
    
    % Find all the lines, and markers
    LineH = findobj(hFig, 'type', 'line', '-or', 'type', 'errorbar');

    if(~isempty(LineH))
        for i=1:length(LineH) % Iterate over all lines in the plot
            % Decide what color for the marker edges
            this_line_color = get(LineH(i),'color');
            if strcmp(marker_outline, 'black')
                marker_outline_color = 'black';
            elseif strcmp(marker_outline, 'matching')
                marker_outline_color = this_line_color;
            else
                marker_outline_color = 'black';
            end

            % If the LineWidth has not been customized, then change it
            if (get(LineH(i), 'LineWidth') <= 1.0)
                set(LineH(i), 'LineWidth', marker_line_width)
            end
            % Change lines and markers if they exist on the plot
            set(LineH(i),   'MarkerSize', marker_size, ...
                'MarkerEdgeColor', marker_outline_color, ...
                'MarkerFaceColor', this_line_color);
        end
    end

    % Find and change the error bars
    LineH = findobj(hFig, 'type', 'errorbar');
    if(~isempty(LineH))
        for i=1:length(LineH) % Iterate over all lines in the plot
            LineH(i).CapSize=error_bar_cap_size;
%             LineH(i).Color = [0 0 0]; % Set all error bars to black

        end
    end

    % Find the legend and subplot title, and if there is one, change it  
    h = get(hFig,'children');
    for k = 1:length(h)
        if strcmpi(get(h(k),'Tag'),'legend')
            set(h(k), 'FontSize', legend_font_size, 'location', 'best');
        end
        if strcmpi(get(h(k),'Type'),'subplottext')
            set(h(k), 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
        end
    end

end