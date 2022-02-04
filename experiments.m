clear;clc;close all

p = parameters();
b = var_bounds(p);

% X = [ 20 10 30;     % outer diameter of float	
%       .3 .1 .5;     % inner diameter ratio of float	
%       30 15 45;     % outer diameter of reaction plate	
%       1 2 3;        % material	
%       10 1 20;      % Number of WECs in array	
%       10 5 50       % D_int	
%       2*pi/7 2*pi/8 2*pi/9
%         ];	

n = 10;
% X = [ b.D_sft_nom, linspace(b.D_sft_min, b.D_sft_max, n);
%       b.D_i_ratio_nom, linspace(b.D_i_ratio_min, b.D_i_ratio_max, n);
%       b.D_or_nom, linspace(b.D_or_min, b.D_or_max, n);
%       1, ones([1 n]);
%       b.N_WEC_nom, linspace(b.N_WEC_min, b.N_WEC_max, n);
%       b.D_int_nom, linspace(b.D_int_min, b.D_int_max, n);
%       b.w_n_nom, linspace(b.w_n_min, b.w_n_max, n) ];

ratios = logspace(-1,1,n);
ratios = [1, ratios(ratios~=1)];
X =  [b.D_sft_nom, b.D_i_ratio_nom, b.D_or_nom, 1, b.N_WEC_nom, b.D_int_nom, b.w_n_nom]' * ratios;
X(4,:) = ones(1,n+1);

X_nom = X(:,1);	
design_size = size(X);	
var_num = design_size(1);	
var_depth = design_size(2);	
LCOE = X*inf;
P_var = X*inf;
opt_idx = zeros(var_num,1);	
recommended = zeros(var_num,2);	
number_runs = var_depth + (var_num-1)*(var_depth-1);	
failed = cell(number_runs,1);	
power = zeros(number_runs,1);	
FOS = zeros(number_runs,1);	
X_ins = zeros(number_runs, var_num);	
design = 0;	
for i = 1:var_num	
    X_in = X_nom;	
    for j = 1:var_depth	
        if i == 1 || j~=1	
            changed_entry = X(i,j);	
            if ~isnan(changed_entry)	
                design = design+1;	
                X_in(i) = changed_entry;	
                X_ins(design,:) = X_in;	
                [LCOE_temp, P_var_temp, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, power(design)] = simulation(X_in,p);	
                FOS(design) = min([FOS1Y,FOS2Y,FOS3Y,FOS_buckling]);
                [feasible, failed{design}] = is_feasible(B, FOS(design), GM, P_var_temp, p);	
                if feasible	
                    LCOE(i,j) = LCOE_temp;	
                    P_var(i,j) = P_var_temp;
                else	
                    LCOE(i,j) = NaN;
                    P_var(i,j) = NaN;
                end	
            end	
        end	
    end	
    [~, opt_idx(i)] = min(LCOE(i,:));	
    recommended(i,:) = [X(i,opt_idx(i)), opt_idx(i)];	
end	
[LCOE_op, ~, ~, B_op, FOS1Y_op,FOS2Y_op,FOS3Y_op,FOS_buckling_op,power_op] = simulation(recommended(:,1),p);	
[op_feasible, CC] = is_feasible(power_op, B_op, FOS1Y_op, power_op, p);	

% create table for display	
var_names = {'D_sft',...    % outer diameter of float (m)	
            'D_i/D_sft',... % inner diameter ratio of float (m)	
            'D_or',...      % outer diameter of reaction plate (m)	
            'M',...         % material (-)	
            'N_WEC',...     % number of WECs in array (-)	
            'D_int',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)
var_names_pretty = {'D_{sft}',...    % outer diameter of float (m)	
            'D_i/D_{sft}',... % inner diameter ratio of float (m)	
            'D_{or}',...      % outer diameter of reaction plate (m)	
            'M',...         % material (-)	
            'N_{WEC}',...     % number of WECs in array (-)	
            'D_{int}',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)
results = array2table(X_ins, 'VariableNames', var_names);	
LCOE = LCOE';
P_var = P_var';
results = addvars(results, round(LCOE(LCOE~=Inf),1), round(power/1e3), round(P_var(P_var~=Inf)), FOS, failed, ...	
    'NewVariableNames', {'LCOE ($/kWh)','Power (kW)','Power Variance (%)','FOS (-)','ConstraintsFailed'});	
disp(results)

% create pareto front
figure
plot(LCOE, P_var, '*')
xlabel('LCOE')
ylabel('P_{var}')
title('Design of Experiments Pareto Front')
legend(var_names_pretty)
improvePlot
