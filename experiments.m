clear;clc;close all

p = parameters();

X = [ 20 10 30;     % outer diameter of float	
      .3 .1 .5;     % inner diameter ratio of float	
      1.5 1.4 1.6;     % outer diameter of reaction plate	
      1 2 3;        % material	
      10 1 20;      % Number of WECs in array	
      10 5 50       % D_int	
      2*pi/7 2*pi/8 2*pi/9
        ];	
    	
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
            'D_or/D_sft',...      % outer diameter ratio of reaction plate (m)	
            'M',...         % material (-)	
            'N_WEC',...     % number of WECs in array (-)	
            'D_int',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)
results = array2table(X_ins, 'VariableNames', var_names);	
LCOE = LCOE';
P_var = P_var';
results = addvars(results, round(LCOE(LCOE~=Inf),1), round(power/1e3), round(P_var(P_var~=Inf)), FOS, failed, ...	
    'NewVariableNames', {'LCOE ($/kWh)','Power (kW)','Power Variance (%)','FOS (-)','ConstraintsFailed'});	
disp(results)