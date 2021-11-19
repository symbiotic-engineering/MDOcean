p = parameters();
b = var_bounds(p);

num_runs = 1e5;
[LCOE,P_var,feasible] = deal(zeros(1,num_runs));

for i=1:num_runs
    x = random_x0(b);

    [LCOE(i), P_var(i), B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM] = simulation(x,p);	
    FOS = min([FOS1Y,FOS2Y,FOS3Y,FOS_buckling]);
    feasible(i) = is_feasible(B, FOS, GM, p);
end
%%
feasible = logical(feasible);
LCOE = LCOE(feasible);
P_var = P_var(feasible);
%%
figure
plot(LCOE,P_var,'*')
hold on
xlabel('LCOE ($/kWh)')
ylabel('Normalized Power Variance (%)')
title('Pareto Front')

[~,idxs] = paretoFront(-1*[LCOE' P_var']);
plot(LCOE(idxs),P_var(idxs),'r*')
xlim([0 10])

plot(fval(:,1),fval(:,2),'ks','MarkerFaceColor','k')

utopia_LCOE = min(LCOE);
utopia_P_var = min(P_var(LCOE<20));
plot(utopia_LCOE,utopia_P_var,'gp','MarkerFaceColor','g','MarkerSize',20)

legend('Dominated Points','Non-Dominated Points','Gradient-Based Pareto Front','Utopia Point')
%%
figure
overallLCOE = [LCOE'; fval(:,1)];
overallPvar = [P_var'; fval(:,2)];
[~,idxo] = paretoFront(-1*[overallLCOE overallPvar]);
plot(overallLCOE(idxo),overallPvar(idxo),'bs','MarkerFaceColor','b')
hold on
plot(min(overallLCOE),min(overallPvar),'gp','MarkerFaceColor','g','MarkerSize',20)

text(2,130, 'Grid-Connected')
text(3,40,  'Microgrids with Storage')
text(6,15, 'Microgrids without Storage')

xlabel('LCOE ($/kWh)')
ylabel('Normalized Power Variance (%)')
title('Combined Pareto Front')
