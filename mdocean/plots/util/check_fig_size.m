function [fig, too_big] = check_fig_size(fig)
    pos = fig.Position(3:4);
    monitor_size = get(groot,'MonitorPositions');
    ratio = pos ./ monitor_size(:,3:4);
    too_big = any(ratio > 1);
    if too_big
        % figure goes offscreen and needs to have its position saved manually to look the same after saving/reopening
        fig.UserData.Position = pos;
    end

end