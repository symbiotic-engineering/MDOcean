p = parameters();
b = var_bounds(p);

num_runs = 1e5;
[LCOE,P_var,feasible] = deal(zeros(1,num_runs));

X = zeros(num_runs,7);
for i=1:num_runs
    x = random_x0(b);
    X(i,:) = x;
    [LCOE(i), P_var(i), B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM] = simulation(x,p);	
    FOS = min([FOS1Y,FOS2Y,FOS3Y,FOS_buckling]);
    feasible(i) = is_feasible(B, FOS, GM, p);
end
%%
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
xlim([0 10])

fval = [11 150]; % this is a dummy value, use fval in workspace from running multiObj.m
plot(fval(:,1),fval(:,2),'ks','MarkerFaceColor','k','HandleVisibility','off')

utopia_LCOE = min(LCOE);
utopia_P_var = min(P_var(LCOE<20));
plot(utopia_LCOE,utopia_P_var,'gp','MarkerFaceColor','g','MarkerSize',20)

legend('Dominated Points','Non-Dominated Points',...'Gradient-Based Pareto Front',
    'Utopia Point')
%%
close all
figure
% overall pareto front
overallLCOE = [LCOE'; fval(:,1)];
overallPvar = [P_var'; fval(:,2)];
[~,idxo] = paretoFront(-1*[overallLCOE overallPvar]);
plot(overallLCOE(idxo),overallPvar(idxo),'bs','MarkerFaceColor','b')
hold on
% utopia point
[minLCOE,idx_best_LCOE] = min(overallLCOE);
[minPvar,idx_best_Pvar] = min(overallPvar);
plot(minLCOE,minPvar,'gp','MarkerFaceColor','g','MarkerSize',20)

% text labels
text(4,30,'RM3 Nominal')
text(minLCOE+.4,minPvar,'Utopia Point')
% text(2,130, 'Grid-Connected')
% text(3,40,  'Microgrids with Storage')
% text(6,15, 'Microgrids without Storage')
text(.3,125,'Solar')

% RM3 nominal reference
x_nom = [b.D_sft_nom, b.D_i_ratio_nom, b.D_or_nom, b.M_nom, b.N_WEC_nom, b.D_int_nom, b.w_n_nom];
[LCOE_nom,P_var_nom] = simulation(x_nom,p);
plot(LCOE_nom,P_var_nom,'rd')

% axis labels
xlabel('LCOE ($/kWh)')
ylabel('Normalized Power Variance (%)')
title('Combined Pareto Front')
xlim([0 10])
improvePlot

% solar reference
plot(0.03, 125,'o','MarkerSize',12,'MarkerEdgeColor',[1, .87, .2],'MarkerFaceColor',[1, .87, .2])
% for the yellow color to work, do not use improvePlot below here

% idenitfy design variables for best designs
x_best_LCOE = X(idx_best_LCOE,:)
x_best_Pvar = X(idx_best_Pvar,:)

% small corner pictures of best geometries
axes('Position',[.25 .7 .15 .15])
box on
visualize_geometry(x_best_LCOE,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])

axes('Position',[.7 .2 .15 .15])
box on
visualize_geometry(x_best_Pvar,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])

axes('Position',[.4 .35 .15 .15])
box on
visualize_geometry(x_nom,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])
