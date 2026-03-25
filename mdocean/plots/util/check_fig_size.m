function [fig, too_big] = check_fig_size(fig)
%CHECK_FIG_SIZE Check if figure exceeds monitor bounds; save position if so
%   If the figure is too large for the monitor, its position is stored in
%   fig.UserData.Position so it can be restored after saving and reopening.
    pos = fig.Position(3:4);
    monitor_size = get(groot,'MonitorPositions');
    ratio = pos ./ monitor_size(:,3:4);
    too_big = any(ratio(:) > 1); % flatten before any() to handle multi-monitor setups
    if too_big
        % figure goes offscreen and needs to have its position saved manually to look the same after saving/reopening
        fig.UserData.Position = pos;
    end
end
