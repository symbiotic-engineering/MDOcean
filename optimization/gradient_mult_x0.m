clear;clc;close all

p = parameters();
b = var_bounds(p);

num_runs = 25;
objs = Inf(num_runs,2);
X_opt = zeros(num_runs,8,2);
flags = zeros(num_runs,2);	
x0s = struct('D_f',[],'D_s_ratio',[],'h_f_ratio',[],'T_s_ratio',[],'F_max',[],'D_int',[],'w_n',[]);

% nominal ICs
[X_opt(1,:,:), objs(1,:), flags(1,:)] = gradient_optim(b.X_nom_struct,p,b);	
x0s(1) = b.X_nom_struct;

% 20 random ICs
for i = 2:num_runs
    [~,x0] = random_x0(b);
    x0s(i) = x0;
    [X_opt(i,:,:), objs(i,:), flags(i,:)] = gradient_optim(x0,p,b);	
end

%% create table for display	
results = struct2table(x0s);
scale = repmat([100 1],num_runs,1);
results = addvars(results, objs.*scale, flags,  ...	
    'NewVariableNames', {'Objs','Flag'});

var_names = {'D_f',...    % outer diameter of float (m)	
            'D_s_ratio',...
            'h_f_ratio',...
            'T_s_ratio',...         
            'F_max',...    
            'D_int',...     % internal damping of controller (Ns/m)	
            'w_n',...       % natural frequency (rad/s)
            'M'};         % material (-)
for i=1:length(var_names)
    if i~= 8
        X = X_opt(:,i,:);
    results = addvars(results, X(:,:), 'NewVariableNames', ...
                    [var_names{i} '_opt'], 'After', var_names{i});
    end
end

results = sortrows(results,{'Flag','Objs'},'descend');
results.Variables =  round(results.Variables,1);
disp(results)
