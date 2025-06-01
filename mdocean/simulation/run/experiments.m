
% Runs one-at-a-time design of experiments

function experiments()

p = parameters();
b = var_bounds();

DOE_strategy = 'ratios'; % 'sample' or 'bounds' or 'ratios'
n = 20;
num_DVs = length(b.X_noms);

if strcmp(DOE_strategy,'sample')
    X = [ 20 10 30;     % outer diameter of float	
          .3 .1 .5;     % inner diameter ratio of float	
          1.5 1.2 2;    % outer diameter of reaction plate	
          1 2 3;        % material	
          10 1 20;      % Number of WECs in array	
          10 5 50       % D_int	
          2*pi/7 2*pi/8 2*pi/9
        ];	
elseif strcmp(DOE_strategy,'bounds')    
    X = zeros(num_DVs,1+n);
    for i = 1:num_DVs
        X(i,1) = b.X_noms(i);
        X(i,2:end) = linspace(b.X_mins(i),b.X_maxs(i),n);
    end
    ratios = X./X(:,1);

elseif strcmp(DOE_strategy,'ratios')
    ratios = logspace(log10(1/3),log10(3),n);
    ratios = [1, ratios(ratios~=1)];
    X =  [b.X_noms] * ratios;
end


X_nom = X(:,1);	
design_size = size(X);	
num_vals_per_DV = design_size(2);
num_vals_swept = num_vals_per_DV - 1; % don't sweep 1st value of each variable (corresponding to ratio 1)

% initialize variables
LCOE = X*inf;
cost = X*inf;
power = X*inf;
opt_idx = zeros(num_DVs,1);	
recommended = zeros(num_DVs,2);	
number_runs = 1 + num_DVs * num_vals_swept; % nominal run plus all sweeps
failed = cell(number_runs,1);	
X_ins = zeros(number_runs, num_DVs);	
design = 0;

% run design of experiments
for i = 1:num_DVs
    X_in = X_nom;	
    for j = 1:num_vals_per_DV
        if i == 1 || j~=1	% prevent rerunning nominal multiple times
            changed_entry = X(i,j);	
            if ~isnan(changed_entry)	
                design = design+1;	
                X_in(i) = changed_entry;	
                X_ins(design,:) = X_in;	

                X_vec = [X_in;b.M_nom];
                [~, ~, failed_lin, feasible_lin] = is_feasible(0, X_vec, p, b);
                
                if feasible_lin
                    [LCOE_temp, cost_temp, P_matrix, g, val] = simulation(X_vec,p);
    
                    % only add to results if relevant nonlin constraints are feasible
                    idx_ignore = false(1,length(b.constraint_names));
                    ignore = {'irrelevant_max_force','LCOE_max','linear_theory','prevent_slamming'};
                    idx_ignore(contains(b.constraint_names,ignore)) = true;
                    [feasible, ~, which_failed] = is_feasible(g, X_in, p, b, idx_ignore);
                else
                    feasible = feasible_lin;
                    which_failed = failed_lin;
                end
                if feasible	
                    LCOE(i,j) = LCOE_temp;	
                    cost(i,j) = cost_temp;
                    power(i,j) = val.power_avg;
                else	
                    LCOE(i,j) = NaN;
                    cost(i,j) = NaN;
                end	
                failed{design} = which_failed;
            end	
        end	
    end	
    [~, opt_idx(i)] = min(LCOE(i,:));	
    recommended(i,:) = [X(i,opt_idx(i)), opt_idx(i)];	
end	

% create table for display	
results = array2table(X_ins, 'VariableNames', b.var_names(1:end-1));	
LCOE = LCOE';
cost = cost';
power = power';
results = addvars(results, round(LCOE(LCOE~=Inf),2), round(cost(cost~=Inf),1), failed, ...
                'NewVariableNames', {'LCOE ($/kWh)','c_v (%)','Failed Constraints'});	
disp(results)

% plot pareto curve for comparison, if pareto results exist
d=dir("**/pareto_search_results*");
if ~isempty(d)
    pareto_curve_heuristics()
    figure(3)
    plot(power/1e3, cost, '*--')
    
    title('Design of Experiments Pareto Front')
    l = legend(b.var_names_pretty);
    improvePlot
    l.Location = 'bestoutside';
end
%% sensitivities plot
[ratios_sorted,idx] = sort(ratios);
LCOE(1,:) = LCOE(1,1); % fill in nominal LCOE results for each DV where it wasn't repeatedly tested
cost(1,:) = cost(1,1);

figure
t = tiledlayout(2,1);
t.TileSpacing = 'compact';

% LCOE subplot
ax1 = nexttile(1);
cols = {'r:','r--','r-','r-.','r.',...       % bulk dims
        'b:','b--',...                       % PTO
        'g:','g--','g-','g-.','g.'};  % structural
yline(LCOE(1,1),'LineWidth',2,'Color','k','HandleVisibility','off')
hold on
for i = 1:size(cols,2)
    temp_LCOE = LCOE(idx,:).';
    plot(ratios_sorted,temp_LCOE(i,:),cols{i})
    hold on
end

ylab1 = ylabel('LCOE ($/kWh)');
x_range = [1/3 3];
axis(ax1,[x_range .6 1.1])
l = legend(b.var_names_pretty{1:end-1});
l.Location = 'northeastoutside';
grid on
hold off

% cost subplot
ax2 = nexttile(2);
yline(cost(1,1),'LineWidth',2,'Color','k')
hold on
for i=1:size(cols,2)
    temp_cost = cost(idx,:).';
    plot(ratios_sorted,temp_cost(i,:),cols{i})
end

ylab2=ylabel('Structural & PTO Cost ($M)');
grid on

% shared plot
title(t,'Design of Experiments Results','FontWeight','bold','FontSize',20)
grid on
linkaxes([ax1,ax2],'x');
xlabel('Design Variable Ratio (-)')
xticklabels(ax1,{})
xticks(ax2,xticks(ax1))
improvePlot
ylab1.FontSize=16.5;
ylab2.FontSize=16.5;
xlim(x_range)
fig = gcf();
fig.Position(3:4) = [600  666]; % make taller

end