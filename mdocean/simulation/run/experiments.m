
% Runs one-at-a-time design of experiments

clear;clc;close all

p = parameters();
b = var_bounds();

DOE_strategy = 'ratios'; % 'sample' or 'bounds' or 'ratios'
n = 20;

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
    X = [ b.D_f_nom, linspace(b.D_f_min, b.D_f_max, n);
          b.D_s_ratio_nom, linspace(b.D_s_ratio_min, b.D_s_ratio_max, n);
          b.h_f_ratio_nom, linspace(b.h_f_ratio_min, b.h_f_ratio_max, n);
          b.T_s_ratio_nom, linspace(b.T_s_ratio_min, b.T_s_ratio_max, n);
          b.F_max_nom, linspace(b.F_max_min, b.F_max_max, n);
          b.B_p_nom, linspace(b.B_p_min, b.B_p_max, n);
          b.w_n_nom, linspace(b.w_n_min, b.w_n_max, n) ];
    ratios = X./X(:,1);
elseif strcmp(DOE_strategy,'ratios')
    ratios = logspace(log10(1/3),log10(3),n);
    ratios = [1, ratios(ratios~=1)];
    X =  [b.D_f_nom b.D_s_ratio_nom b.h_f_ratio_nom b.T_s_ratio_nom b.F_max_nom b.B_p_nom b.w_n_nom b.M_nom]' * ratios;
end

X_nom = X(:,1);	
design_size = size(X);	
num_vars = design_size(1);	
num_vals_per_var = design_size(2);
num_vals_swept = num_vals_per_var - 1; % don't sweep 1st value of each variable (corresponding to ratio 1)

% don't sweep material
idx_material = 8;
X(idx_material,:) = b.M_nom * ones(1,n+1);
num_vars_swept = num_vars - 1;

% initialize variables
LCOE = X*inf;
P_var = X*inf;
opt_idx = zeros(num_vars,1);	
recommended = zeros(num_vars,2);	
number_runs = 1 + num_vars_swept * num_vals_swept; % nominal run plus all sweeps
failed = cell(number_runs,1);	
power = zeros(number_runs,1);	
FOS = zeros(number_runs,1);	
X_ins = zeros(number_runs, num_vars);	
design = 0;

% run design of experiments
for i = 1:num_vars_swept
    X_in = X_nom;	
    for j = 1:num_vals_per_var
        if i == 1 || j~=1	% prevent rerunning nominal multiple times
            changed_entry = X(i,j);	
            if ~isnan(changed_entry)	
                design = design+1;	
                X_in(i) = changed_entry;	
                X_ins(design,:) = X_in;	
                [LCOE_temp, P_var_temp, P_matrix, g, val] = simulation(X_in,p);

                % only add to results if first 12 constraints are feasible
                g(13) = 0;
                g(14) = 0;
                [feasible, which_failed] = is_feasible(g, X_in, p, b);
                if feasible	
                    LCOE(i,j) = LCOE_temp;	
                    P_var(i,j) = P_var_temp;
                else	
                    LCOE(i,j) = NaN;
                    P_var(i,j) = NaN;
                end	
                failed{design} = which_failed;
            end	
        end	
    end	
    [~, opt_idx(i)] = min(LCOE(i,:));	
    recommended(i,:) = [X(i,opt_idx(i)), opt_idx(i)];	
end	

% create table for display	
results = array2table(X_ins, 'VariableNames', b.var_names);	
LCOE = LCOE';
P_var = P_var';
results = addvars(results, round(LCOE(LCOE~=Inf),2), round(P_var(P_var~=Inf),1), failed, ...
                'NewVariableNames', {'LCOE ($/kWh)','c_v (%)','Failed Constraints'});	
disp(results)

% plot pareto curve for comparison
pareto_curve_heuristics()
figure(2)
plot(LCOE, P_var, '*--')
xlabel('LCOE')
ylabel('P_{var}')
title('Design of Experiments Pareto Front')
l = legend(b.var_names_pretty);
improvePlot
l.Location = 'bestoutside';

%% sensitivities plot
[ratios_sorted,idx] = sort(ratios);
LCOE(1,:) = LCOE(1,1); % fill in nominal LCOE results for each DV where it wasn't repeatedly tested
P_var(1,:) = P_var(1,1);

figure
t = tiledlayout(2,1);
t.TileSpacing = 'compact';

ax1 = nexttile(1);
plot(ratios_sorted,LCOE(idx,:).')
ylabel('LCOE ($/kWh)')
grid on
ax2 = nexttile(2);
plot(ratios_sorted,P_var(idx,:).')
ylabel('Power c_v (%)')
grid on

title(t,'Design of Experiments Results','FontWeight','bold','FontSize',20)
l = legend(b.var_names_pretty);
l.Location = 'bestoutside';
grid on
linkaxes([ax1,ax2],'x');
xlabel('Design Variable Ratio (-)')
xticklabels(ax1,{})
xticks(ax2,xticks(ax1))
improvePlot

