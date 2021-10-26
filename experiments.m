clear;clc;close all

p = parameters();

X = [ 20 10 30;     % outer diameter of float
      .3 .1 .5;     % inner diameter ratio of float
      30 15 45;     % outer diameter of reaction plate
      1 2 3;        % material
      10 1 20;      % Number of WECs in array
      1e7 5e6 5e7   % D_int
        ];
    
X_nom = X(:,1);
design_size = size(X);
var_num = design_size(1);
var_depth = design_size(2);
LCOE = X*inf;
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
                [LCOE_temp, ~, ~, B, FOS(design), power(design)] = simulation(X_in,p);
                [feasible, failed{design}] = is_feasible(power(design), B, FOS(design), p);
                if feasible
                    LCOE(i,j) = LCOE_temp;
                else
                    LCOE(i,j) = NaN;
                end
            end
        end
    end
    [~, opt_idx(i)] = min(LCOE(i,:));
    recommended(i,:) = [X(i,opt_idx(i)), opt_idx(i)];
end
[LCOE_op, ~, ~, B_op, FOS_op, power_op] = simulation(recommended(:,1),p);
[op_feasible, CC] = is_feasible(power_op, B_op, FOS_op, p);

% create table for display
var_names = {'D_sft',...    % outer diameter of float (m)
            'D_i/D_sft',... % inner diameter ratio of float (m)
            'D_or',...      % outer diameter of reaction plate (m)
            'M',...         % material (-)
            'N_WEC',...     % number of WECs in array (-)
            'D_int'};       % internal damping of controller (Ns/m)
results = array2table(X_ins, 'VariableNames', var_names);
LCOE = LCOE';
results = addvars(results, round(LCOE(LCOE~=Inf),1), round(power/1e3), FOS, failed, ...
    'NewVariableNames', {'LCOE ($/kWh)','Power (kW)','FOS (-)','ConstraintsFailed'});
disp(results)
