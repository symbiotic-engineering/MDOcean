p = parameters();
b = var_bounds(p);

num_runs = 1e5;
[LCOE,P_var,feasible] = deal(zeros(1,num_runs));

X = zeros(num_runs,7);
for i=1:num_runs
    x = random_x0(b);
    X(i,:) = x;
    [LCOE(i), P_var(i), B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM, P_elec] = simulation(x,p);	
    FOS = min([FOS1Y,FOS2Y,FOS3Y,FOS_buckling]);
    feasible(i) = is_feasible(B, FOS, GM, P_elec, p);
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
LCOE_max = 1;
xlim([0 LCOE_max])

fval = [11 150]; % this is a dummy value, use fval in workspace from running multiObj.m
plot(fval(:,1),fval(:,2),'ks','MarkerFaceColor','k','HandleVisibility','off')

utopia_LCOE = min(LCOE);
utopia_P_var = min(P_var(LCOE<LCOE_max));
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
filteredOverallPvar = overallPvar(overallLCOE<LCOE_max);
filteredOverallLCOE = overallLCOE(overallLCOE<LCOE_max);
[minPvar,idx_best_Pvar] = min(filteredOverallPvar);
plot(minLCOE,minPvar,'gp','MarkerFaceColor','g','MarkerSize',20)

% balanced design
[~,idx_balanced] = min(abs(overallPvar-100));
LCOE_balanced = overallLCOE(idx_balanced);
P_var_balanced = overallPvar(idx_balanced);

% RM3 nominal reference
x_nom = [b.D_sft_nom, b.D_i_ratio_nom, b.D_or_nom, b.M_nom, b.N_WEC_nom, b.D_int_nom, b.w_n_nom]
[LCOE_nom,P_var_nom] = simulation(x_nom,p);
plot(LCOE_nom,P_var_nom,'rd')

% axis labels
xlabel('LCOE ($/kWh)')
ylabel('Power Variation (%)')
title('Pareto Front')
xlim([0 LCOE_max])
ylim([50 275])
improvePlot

% solar reference
LCOE_solar = 0.03;
P_var_solar = 125;
plot(LCOE_solar, P_var_solar,'o','MarkerSize',12,'MarkerEdgeColor',[1, .87, .2],'MarkerFaceColor',[1, .87, .2])
% for the yellow color to work, do not use improvePlot below here

% text labels
text(LCOE_nom+.03,P_var_nom,'RM3 Nominal')
text(minLCOE+.03,minPvar,'Utopia Point')
% text(2,130, 'Grid-Connected')
% text(3,40,  'Microgrids with Storage')
% text(6,15, 'Microgrids without Storage')
text(LCOE_solar+.03,P_var_solar,'Solar')
text(overallLCOE(idx_best_LCOE)+.03,overallPvar(idx_best_LCOE),'Cheapest')
text(filteredOverallLCOE(idx_best_Pvar)-.2,filteredOverallPvar(idx_best_Pvar),'Least Variable')
text(LCOE_balanced+.04,P_var_balanced+5,'Balanced Design')

% idenitfy design variables for best designs
x_best_LCOE = X(idx_best_LCOE,:)
x_best_Pvar = X(idx_best_Pvar,:)
x_balanced = X(idx_balanced,:)

% small corner pictures of best geometries
% upper left
axes('Position',[.15 .65 .15 .2])
box on
visualize_geometry(x_best_LCOE,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])

% lower right
axes('Position',[.7 .25 .15 .2])
box on
visualize_geometry(x_best_Pvar,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])

% balanced
axes('Position',[.38 .33 .15 .2])
box on
visualize_geometry(x_balanced,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])

% RM3
axes('Position',[.55 .63 .15 .2])
box on
visualize_geometry(x_nom,p,true);
set(gca,'XTickLabel',[],'YTickLabel',[])