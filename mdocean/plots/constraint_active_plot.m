function [idx] = constraint_active_plot(residuals,fval,tol,b,reversed)

    if nargin<5
        reversed = false;
    end

    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    nlcon_active = abs(residuals.ineqnonlin) < tol;
    lincon_active = abs(residuals.ineqlin) < tol;

    if reversed
        dir = 'descend';
    else
        dir = 'ascend';
    end
    [~,idx] = sort(fval(:,1),dir); % order by increasing LCOE (decreasing if reversed)

    figure
    tiledlayout(2,2); 
    
    nexttile
    spy(lb_active(idx,:)');
    title('Lower Bound Active')
    grid on
    set(gca,'ytick',1:size(lb_active,2),'yticklabel',b.var_names_pretty(1:end-1))
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])

    nexttile
    spy(ub_active(idx,:)')
    title('Upper Bound Active')
    grid on
    set(gca,'ytick',1:size(ub_active,2),'yticklabel',b.var_names_pretty(1:end-1))
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])

    nexttile
    spy(nlcon_active(idx,:)')
    title('Nonlinear Constraint Active')
    grid on
    set(gca,'ytick',1:size(nlcon_active,2),'yticklabel',b.constraint_names_pretty)
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])

    nexttile
    spy(lincon_active(idx,:)')
    title('Linear Constraint Active')
    set(gca,'ytick',1:size(lincon_active,2),'yticklabel',b.lin_constraint_names_pretty)
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1])
    grid on

    %improvePlot
end