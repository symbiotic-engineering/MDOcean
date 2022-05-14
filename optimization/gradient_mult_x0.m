clear;clc;close all

p = parameters();
b = var_bounds(p);

num_runs = 20;
LCOE = Inf(num_runs,1);
X_opt = zeros(num_runs,7);
flag = zeros(num_runs,1);	
x0s = struct('D_sft',[],'D_i_ratio',[],'D_or_ratio',[],'N_WEC',[],'D_int',[],'w_n',[]);

for i = 1:num_runs
    [~,x0] = random_x0(b);
    x0s(i) = x0;
    [X_opt(i,:), LCOE(i), flag(i)] = gradient_optim(x0,p,b);	
end

% create table for display	
var_names = {'D_sft',...    % outer diameter of float (m)	
            'D_i/D_sft',... % inner diameter ratio of float (m)	
            'D_or/D_sft',...      % outer diameter of reaction plate (m)	
            'M',...         % material (-)	
            'N_WEC',...     % number of WECs in array (-)	
            'D_int',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)
        %%
results = struct2table(x0s);
results = addvars(results, round(LCOE,1), flag,  ...	
    'NewVariableNames', {'LCOE ($/kWh)','Flag'});

var_names = {'D_sft',...    % outer diameter of float (m)	
            'D_i_ratio',... % inner diameter ratio of float (m)	
            'D_or_ratio',...      % outer diameter of reaction plate (m)	
            'M',...         % material (-)	
            'N_WEC',...     % number of WECs in array (-)	
            'D_int',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)
for i=1:length(var_names)
    if i~= 4
    results = addvars(results, X_opt(:,i), 'NewVariableNames', ...
                    [var_names{i} '_opt'], 'After', var_names{i});
    end
end

results = sortrows(results,'Flag','descend');
results.Variables =  round(results.Variables,1);
disp(results)

