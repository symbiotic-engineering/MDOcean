function [] = delta_x(X,grad,hess,J,p,b,which_obj)
%DELTA_X Use grad to estimate delta J for given delta x

% capital X has all DVs, lowercase x has only continuous DVs

x = X(1:end-1);
dx = x/10;

% don't use the hessian contribution because it's the hess of the
% Lagrangian, not of the objective
dJ_up   =  grad .* dx + 0*1/2 * dx.^2 .* diag(hess);
dJ_down = -grad .* dx + 0*1/2 * dx.^2 .* diag(hess);

J_up   = zeros(1,length(x));
J_down = zeros(1,length(x));
for i=1:length(x)
    dX = zeros(size(X));
    dX(i) = dx(i);
    [out_up(1),  out_up(2)]   = simulation(X+dX,p);
    J_up(i) = out_up(which_obj);
    [out_down(1),out_down(2)] = simulation(X-dX,p);
    J_down(i) = out_down(which_obj);
end

dJ_up_actual = J_up - J;
dJ_down_actual = J_down - J;

f = figure;
dJ = [dJ_up';dJ_up_actual;dJ_down';dJ_down_actual];
var = strcat('$',b.var_names_pretty(1:end-1),'$');
var_cats = categorical(var,var,var);
bar(var_cats, dJ/J)
title('$\frac{\Delta J}{J}$ for $\Delta x = \pm x^*/10$','Interpreter','latex')
legend('$+\Delta x$: Gradient Estimate','$+\Delta x$: Actual',...
       '$-\Delta x$: Gradient Estimate','$-\Delta x$: Actual','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$\Delta J/J$','Interpreter','latex')
improvePlot
set(gca(),'TickLabelInterpreter','latex','YGrid','on')
ylim([-.09 .09])
f.Position(3:4) = [705 475];

end

