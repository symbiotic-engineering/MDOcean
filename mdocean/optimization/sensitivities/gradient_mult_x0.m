function [results] = gradient_mult_x0(filename_uuid)

p = parameters();
b = var_bounds();

if nargin>0
    b.filename_uuid = filename_uuid;
end

num_runs = 3;
num_DVs = length(b.var_names);
num_objs = 2;
objs = Inf(num_runs,num_objs);
X_opt = zeros(num_runs,num_DVs,num_objs);
flags = zeros(num_runs,num_objs);	

% nominal ICs
[X_opt(1,:,:), objs(1,:), flags(1,:)] = gradient_optim(b.X_start_struct,p,b);	
x0s(1) = b.X_start_struct;

% 20 random ICs
for i = 2:num_runs
    [x0_vec,x0] = random_x0(b);
    [~, ~, feasible_lin] = is_feasible(0, x0_vec, p, b);
    if feasible_lin
        x0s(i) = x0;
        [X_opt(i,:,:), objs(i,:), flags(i,:)] = gradient_optim(x0,p,b);	
    end
end

%% create table for display	

results = struct2table(x0s);
cents_per_dollar = 100;
scale = repmat([cents_per_dollar 1],num_runs,1); % scale LCOE units
results = addvars(results, objs.*scale, flags,  ...	
    'NewVariableNames', {'Objs','Flag'});

for i=1:length(b.var_names)-1
    X = X_opt(:,i,:);
    results = addvars(results, X(:,:), 'NewVariableNames', ...
                    [b.var_names{i} '_opt'], 'After', b.var_names{i});
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

end