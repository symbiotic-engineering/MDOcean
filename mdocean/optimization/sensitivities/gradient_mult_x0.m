clear;clc;close all

p = parameters();
b = var_bounds();

num_runs = 100;
objs = Inf(num_runs,2);
X_opt = zeros(num_runs,8,2);
flags = zeros(num_runs,2);	
x0s = struct('D_f',[],'D_s_ratio',[],'h_f_ratio',[],'T_s_ratio',[],'F_max',[],'D_int',[],'w_n',[]);

% nominal ICs
[X_opt(1,:,:), objs(1,:), flags(1,:)] = gradient_optim(b.X_start_struct,p,b);	
x0s(1) = b.X_start_struct;

% 20 random ICs
for i = 2:num_runs
    [~,x0] = random_x0(b);
    x0s(i) = x0;
    [X_opt(i,:,:), objs(i,:), flags(i,:)] = gradient_optim(x0,p,b);	
end

%% create table for display	
[x0s(1:num_runs).B_p] = deal(x0s.D_int);
results = struct2table(x0s);
scale = repmat([100 1],num_runs,1);
results = addvars(results, objs.*scale, flags,  ...	
    'NewVariableNames', {'Objs','Flag'});

for i=1:length(b.var_names)
    if i~= 8
        X = X_opt(:,i,:);
    results = addvars(results, X(:,:), 'NewVariableNames', ...
                    [b.var_names{i} '_opt'], 'After', b.var_names{i});
    end
end

results = sortrows(results,{'Flag','Objs'},'descend');
results.Variables =  round(results.Variables,1);
disp(results)

%% summary statistics
converged = flags > 0;
kkt = flags == 1;
tol = 0.001;
optimal = abs(objs - min(objs)) < tol; % this assumes that solutions converging 
% to the same objective have the same optimal DVs, which is not guaranteed 
% but seems correct in this case by inspection of the table above

optimal_and_converged = converged & optimal;
optimal_and_kkt = kkt & optimal;

percent_converged = sum(converged) / num_runs
percent_kkt = sum(kkt) / num_runs
percent_optimal = sum(optimal) / num_runs
percent_optimal_given_converged = sum(optimal_and_converged) ./ sum(converged)
percent_optimal_given_kkt = sum(optimal_and_kkt) ./ sum(kkt)
