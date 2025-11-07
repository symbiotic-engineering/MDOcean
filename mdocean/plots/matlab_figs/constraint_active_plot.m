function [idx,f] = constraint_active_plot(residuals,fval,tol,b,reversed)

    if nargin<5
        reversed = false;
    end

    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    nlcon_active = abs(residuals.ineqnonlin) < tol;
    lincon_active = abs(residuals.ineqlin) < tol;

    % merge sea state slamming constraints
    idx_slamming = contains(b.constraint_names,'slamming');
    idx_slamming_first = strcmp(b.constraint_names,'prevent_slamming1');
    idx_slamming_after = idx_slamming & ~idx_slamming_first;
    nlcon_active(:,idx_slamming_first) = any(nlcon_active(:,idx_slamming),2);
    nlcon_active(:,idx_slamming_after) = [];
    constraint_names_mod = b.constraint_names_pretty(~idx_slamming_after);
    constraint_names_mod(idx_slamming_first) = {'Prevent Slamming'};

    if reversed
        dir = 'descend';
    else
        dir = 'ascend';
    end
    [~,idx] = sort(fval(:,1),dir); % order by increasing LCOE (decreasing if reversed)

    f = figure('Color','w');
    tiledlayout(2,3,'TileSpacing','tight','Padding','tight');

    nexttile
    spy(lb_active(idx,:)','v');
    hold on
    spy(ub_active(idx,:)','^r')
    title('Design Variable Bound Active','FontSize',14)
    grid on
    x_string = strcat('x_{', num2str((1:length(b.var_names)-1).') , "} ", b.var_names_pretty(1:end-1).' );
    set(gca,'ytick',1:size(lb_active,2),'yticklabel',x_string)
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])
    xlabel('Number Along Pareto Front','FontSize',14)
    legend('Lower bound','Upper bound')

    nexttile(4)
    spy(lincon_active(idx,:)')
    title('Linear Constraint Active','FontSize',14)
    lin_string = strcat('g_{L,', num2str((1:length(b.lin_constraint_names)).') , "} ", b.lin_constraint_names_pretty.' );
    set(gca,'ytick',1:size(lincon_active,2),'yticklabel',lin_string)
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])
    grid on
    xlabel('Number Along Pareto Front','FontSize',14)

    nexttile(2,[2 2])
    spy(nlcon_active(idx,:)')
    title('Nonlinear Constraint Active','FontSize',14)
    grid on
    nl_string = strcat('g_{NL,', num2str((1:length(constraint_names_mod)).') , "} ", constraint_names_mod.' );
    set(gca,'ytick',1:size(nlcon_active,2),'yticklabel',nl_string)
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])
    xlabel('Number Along Pareto Front','FontSize',14)

    f.Position = [2 128 1532  626];
    %improvePlot
end
