function fig = lagrange_multiplier_bar_chart(b,lambda)
%LAGRANGE_MULTIPLIER_BAR_CHART 

idx_lb         = lambda.lower ~= 0;
idx_ub         = lambda.upper ~= 0;
idx_ineqlin    = lambda.ineqlin ~= 0;
%idx_eqlin      = lambda.eqlin ~= 0;
idx_ineqnonlin = lambda.ineqnonlin ~= 0;
%idx_eqnonlin   = lambda.eqnonlin ~= 0;

x_lb = strcat(b.var_names_pretty(idx_lb), ' Lower Bound');
x_ub = strcat(b.var_names_pretty(idx_ub), ' Upper Bound');
x_ineqlin = b.lin_constraint_names_pretty(idx_ineqlin);
x_ineqnonlin = b.constraint_names_pretty(idx_ineqnonlin);

X = [x_lb x_ub x_ineqlin x_ineqnonlin];
X = reordercats(categorical(X),X);

y_lb = lambda.lower(idx_lb);
y_ub = lambda.upper(idx_ub);
y_ineqlin = lambda.ineqlin(idx_ineqlin);
y_ineqnonlin = lambda.ineqnonlin(idx_ineqnonlin);
Y = [y_lb; y_ub; y_ineqlin; y_ineqnonlin];

fig = figure;
bar(X,Y)
title('Lagrange Multipliers \lambda=-\partial J/\partial g')
improvePlot

end

