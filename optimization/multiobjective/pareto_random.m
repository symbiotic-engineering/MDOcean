clear;clc
p = parameters();
b = var_bounds(p);
num_runs = 1e5;
[LCOE,P_var,feasible] = deal(zeros(1,num_runs));

X = zeros(num_runs,8);
for i=1:num_runs
    xx = random_x0(b);
    X(i,:) = xx;
    [LCOE(i), P_var(i), B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, D_d, P_matrix, g] = simulation(xx, p);	
    FOS = min([FOS1Y,FOS2Y,FOS3Y,FOS_buckling]);
    feasible(i) = is_feasible(B, FOS, GM, P_elec, D_d, g(16), g(17), g(18), p);
end

feasible = logical(feasible);
LCOE(~feasible) = Inf;
P_var(~feasible) = Inf;
%%
figure
plot(LCOE,P_var,'*')
hold on
xlabel('LCOE ($/kWh)')
ylabel('Normalized Power Variance (%)')
title('Pareto Front')

[~,idxs] = paretoFront(-1*[LCOE' P_var']);
plot(LCOE(idxs),P_var(idxs),'r*')
LCOE_max = p.LCOE_max;

utopia_LCOE = min(LCOE);
utopia_P_var = min(P_var(LCOE<LCOE_max));
plot(utopia_LCOE,utopia_P_var,'gp','MarkerFaceColor','g','MarkerSize',20)

legend('Dominated Points','Non-Dominated Points','Utopia Point')