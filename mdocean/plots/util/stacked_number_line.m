
% uncomment to test
% rng('default')
% rand_data = rand(5,4);
% low = [1 2 3 4 5];
% up = [6 7 8 9 10];

% data = zeros(5,4);
% for i=1:5
%     range = up(i) - low(i);
%     data(i,:) = rand_data(i,:) * range + low(i);
% end
% figure
% stacked_number_line(data, low, up, {'r','b','g','k'},{'o','x','^','v'},{'a','b','c','d'},{'1','2','3','4','5'})

function stacked_number_line(data, lower_lim, upper_lim, colors, markers, legend_cell, ylabel_cell)
%STACKED_NUMBER_LINE Plot each row of a matrix as a number line.
%
% stacked_number_line(data)
% stacked_number_line(data, lower_lim, upper_lim)
% stacked_number_line(data, lower_lim, upper_lim, colors, markers, legend_cell, ylabel_cell)
%
% Inputs
% -------
% data : NxM numeric matrix
%     Each row is plotted as markers at (x,y) = (data(i,:), 0).
%
% lower_lim : Nx1 or 1xN numeric vector (optional)
%     Lower x-axis limit for each row.
%
% upper_lim : Nx1 or 1xN numeric vector (optional)
%     Upper x-axis limit for each row.
%
% colors : 1xM cell array of color specifiers (optional)
%     Color for each marker series.
%
% markers : 1xM cell array of marker specifiers (optional)
%     Marker type for each marker series.
%
% legend_cell : 1xM cell array of character vectors or strings (optional)
%     Legend labels for marker series.
%
% ylabel_cell : Nx1 or 1xN cell array of character vectors or strings (optional)
%     Y-axis labels for each row.
%
% Notes
% -----
% - Only the x-axis is shown.
% - Markers are plotted without connecting lines.
% - If limits are omitted, MATLAB chooses limits automatically.

    arguments
        data (:,:) double
        lower_lim (:,1) double = []
        upper_lim (:,1) double = []
        colors = {}
        markers = {}
        legend_cell = {}
        ylabel_cell = {}
    end

    [N,M] = size(data);

    % Accept row or column vectors
    if ~isempty(lower_lim)
        lower_lim = lower_lim(:);
        assert(numel(lower_lim)==N, ...
            'lower_lim must have length N.');
    end

    if ~isempty(upper_lim)
        upper_lim = upper_lim(:);
        assert(numel(upper_lim)==N, ...
            'upper_lim must have length N.');
    end

    t = tiledlayout(N,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for i = 1:N
        ax = nexttile(t);

        % Plot markers only
        for j = 1:M
            style = [colors{j} markers{j}];
            plot(ax, data(i,j), 0, style, 'LineStyle','none');
            hold(ax,'on');
        end

        % Hide y-axis entirely
        ax.YTick = [];
        ax.YColor = 'none';
        ax.Box = 'off';

        % Keep only bottom x-axis
        ax.XAxisLocation = 'origin';

        % Small vertical extent so markers are visible
        ylim(ax,[-0.1 0.1]);

        % Remove unnecessary y margins
        ax.YLimMode = 'manual';

        % Set x limits if provided
        if ~isempty(lower_lim) && ~isempty(upper_lim)
            lims = [lower_lim(i), upper_lim(i)];
        elseif ~isempty(lower_lim)
            xl = xlim(ax);
            lims = [lower_lim(i), xl(2)];
        elseif ~isempty(upper_lim)
            xl = xlim(ax);
            lims = [xl(1), upper_lim(i)];
        else
            lims = xlim(ax);
        end
        
        axis padded
        if i==1
            leg = legend(legend_cell);
            leg.Layout.Tile = 'east';
        end
        xlim(ax,lims)
        ax.XAxis.Color = 'w';

        %ax.XAxis.TickLabelColor = 'g'; to debug and ensure axes overlap ok
        oldpos = ax.Position;
        newpos = oldpos;
        newpos(4) = oldpos(4)/2; % half height
        ax2 = axes('position', newpos,'Color','none','XAxisLocation','top');
        ylabel(ax,ylabel_cell{i})
        ax2.YAxis.Color = 'none';
        ax.YLabel.Color = 'k';
        xlim(ax2,lims)
    end
end

