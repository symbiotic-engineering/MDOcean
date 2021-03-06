clear;clc
p = parameters();
b = var_bounds(p);

num_runs = 1e4;
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
xlim([0 LCOE_max])

%[x,fval] = pareto_search();
load("pareto_search_results6.mat")
cols = [1 3 6 5 4 2 7];
X_ps = x(:,cols); % swap indices based on solver generated function
X_ps = [X_ps ones(length(X_ps),1)]; % add eigth column for material 

plot(fval(:,1),fval(:,2),'ks','MarkerFaceColor','k','HandleVisibility','off')

utopia_LCOE = min(LCOE);
utopia_P_var = min(P_var(LCOE<LCOE_max));
plot(utopia_LCOE,utopia_P_var,'gp','MarkerFaceColor','g','MarkerSize',20)

legend('Dominated Points','Non-Dominated Points','Pareto Search Results',...
    'Utopia Point')
%%
simplePareto = true; % toggle between a simple pareto front and one that's 
% annotated with the three recommended designs

%close all
figure
% overall pareto front
overallLCOE = [LCOE'; fval(:,1)];
overallPvar = [P_var'; fval(:,2)];
overallX = [X; X_ps];
[~,idxo] = paretoFront(-1*[overallLCOE overallPvar]);
plot(overallLCOE(idxo),overallPvar(idxo),'bs','MarkerFaceColor','b')
hold on

% utopia point
[minLCOE,idx_best_LCOE] = min(overallLCOE);
[minPvar,idx_best_Pvar] = min(overallPvar);
plot(minLCOE,minPvar,'gp','MarkerFaceColor','g','MarkerSize',20)

% RM3 nominal reference
x_nom = b.X_noms;
x_nom(8) = 1;
[LCOE_nom,P_var_nom] = simulation(x_nom,p);
LCOE_nom = 0.75; % from RM3 report p175 fig 5-33
plot(LCOE_nom,P_var_nom,'rd')

if ~simplePareto
    % balanced design
    [~,idx_balanced] = min(abs(overallPvar-100));
    LCOE_balanced = overallLCOE(idx_balanced);
    P_var_balanced = overallPvar(idx_balanced);
    
    % black squares for 3 ref points
    plot(minLCOE,overallPvar(idx_best_LCOE),'ks')
    plot(overallLCOE(idx_best_Pvar),minPvar,'ks')
    %plot(LCOE_balanced,P_var_balanced,'ks')     
    plot(0.18,100,'ks') % hardcode solution from running gradient_optim with maxLCOE set to 0.18
end

% axis labels
xlabel('LCOE ($/kWh)')
ylabel('Power Variation (%)')
title('Pareto Front')
xlim([0 LCOE_max])
ylim([50 290])
improvePlot

% solar reference
LCOE_solar = 0.03;
P_var_solar = 125;
plot(LCOE_solar, P_var_solar,'o','MarkerSize',12,'MarkerEdgeColor',[1, .87, .2],'MarkerFaceColor',[1, .87, .2])
% for the yellow color to work, do not use improvePlot below here

% text labels
text(LCOE_nom+.03,P_var_nom,'RM3 Nominal')
text(minLCOE+.03,minPvar,'Utopia Point')
text(LCOE_solar,P_var_solar-10,'Solar')
if ~simplePareto
    text(overallLCOE(idx_best_LCOE)+.03,overallPvar(idx_best_LCOE),'Cheapest')
    text(overallLCOE(idx_best_Pvar)-.2,overallPvar(idx_best_Pvar),'Least Variable')
    text(LCOE_balanced-.04,P_var_balanced+5,'Balanced Design')
end
% idenitfy design variables for best designs
x_best_LCOE = overallX(idx_best_LCOE,:)
x_best_Pvar = overallX(idx_best_Pvar,:)
%x_balanced = overallX(idx_balanced,:)

x_balanced = [13.0468, 0.4599, 0.1000, 0.9738, 12.7445, 30.6899, 10.1115, 1.0000]; % hardcode from gradient optim

if ~simplePareto
    % small corner pictures of best geometries
    % upper left
    axes('Position',[.15 .7 .15 .2])
    box on
    visualize_geometry(x_best_LCOE,p,true);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    
    % lower right
    axes('Position',[.7 .2 .15 .2])
    box on
    visualize_geometry(x_best_Pvar,p,true);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    
    % balanced
    axes('Position',[.33 .33 .15 .2])
    box on
    visualize_geometry(x_balanced,p,true);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    
    % RM3
    axes('Position',[.55 .63 .15 .2])
    box on
    visualize_geometry(x_nom,p,true);
    set(gca,'XTickLabel',[],'YTickLabel',[])
end

%% plots for DVs as a fn of percent along the pareto
LCOE_pareto = overallLCOE(idxo);
Pvar_pareto = overallPvar(idxo);
[LCOE_pareto_sorted,idx_sort] = sort(LCOE_pareto(LCOE_pareto<LCOE_max));
Pvar_pareto_sorted = Pvar_pareto(LCOE_pareto<LCOE_max);
Pvar_pareto_sorted = Pvar_pareto_sorted(idx_sort);

pct = linspace(0,100,length(idx_sort));
%pct_angle = 100/(pi/2) * atan((LCOE_pareto_sorted - minLCOE) ./ (Pvar_pareto_sorted - minPvar));
num = (overallPvar(idx_best_LCOE) - Pvar_pareto_sorted) / (overallPvar(idx_best_LCOE) - minPvar);
den = (overallLCOE(idx_best_Pvar) - LCOE_pareto_sorted) / (overallLCOE(idx_best_Pvar) - minLCOE);
pct_angle = 100/(pi/2) * atan( num ./ den);
pct_angle(pct_angle==-100) = 100;

X_pareto = overallX(idxo,:);
X_pareto = X_pareto(LCOE_pareto<LCOE_max,:);
X_pareto_sorted = X_pareto(idx_sort,:);

X_pareto_sorted_scaled = X_pareto_sorted ./ repmat(x_best_LCOE,length(idx_sort),1);

X_pareto_sorted_scaled = X_pareto_sorted_scaled(:,1:7); % get rid of material

windowSize = round(length(idx_sort) * 5/100);
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = zeros(size(X_pareto_sorted_scaled));
for i=1:7
    x = X_pareto_sorted_scaled(:,i); 
    x_padded = [ones(windowSize,1); x];
    yy = filter(b,a,x_padded); % moving average filter
    y(:,i) = yy((windowSize+1):end);
end

var_names_pretty = {'D_f',...    % outer diameter of float (m)	
            'D_s/D_f',...
            'h_f/D_f',...      	
            'T_s/h_s',...      	
            'F_{max}',...     % max force (N)
            'D_{int}',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)


% unfiltered
figure
semilogy(pct_angle,X_pareto_sorted_scaled)
title('Unfiltered Design Heuristics')
xlabel('Percent along the Pareto Curve')
ylabel('Normalized Optimal Design Value')
legend(var_names_pretty,'Location','eastoutside')
ylim([0 15])
improvePlot
grid on
set(gca,'YMinorGrid','on')

% filtered
cols = {'r:','r--','r-','r-.','b:','b--','b-'};
figure
for i=1:7
    semilogy(pct_angle,y(:,i),cols{i})
    hold on
end

% make fake major grid lines (didn't do grid on because then major lines show up for .2 .5 2 5 too)
x_grid = [3 97];
y_grid = [.1 1 10 100 1000];
for i=1:length(y_grid)
    plot(x_grid, y_grid(i)*[1 1],'Color',[.85 .85 .85]);
end

title('Design Heuristics')
%xlabel('Percent along the Pareto Curve')
ylabel('Normalized Optimal Design Value')
ylim([.03 5000])
set(gca,'YTick',y_grid)
improvePlot
legend(var_names_pretty,'Location','eastoutside')
set(gca,'YMinorGrid','on')
set(gca,'XGrid','on')
set(gca, 'Children', flipud(get(gca, 'Children')) ) % put fake gridlines behind real lines

figure
plot(pct_angle,LCOE_pareto_sorted*100,pct_angle,Pvar_pareto_sorted)
grid on
xlabel('Percent along the Pareto Curve')
ylabel('Objective Value')
improvePlot
cent = char(0162);
legend(['LCOE (' cent '/kWh)'],'c_v (%)',Location='northeast')
