function [results] = gradient_mult_x0(filename_uuid)

p = parameters();
b = var_bounds();

if nargin>0
    b.filename_uuid = filename_uuid;
end

num_runs = 1000;
num_DVs = length(b.var_names);
num_objs = 2;
objs = Inf(num_runs,num_objs);
X_opt = zeros(num_runs,num_DVs,num_objs);
flags = zeros(num_runs,num_objs);	

% nominal ICs
[X_opt(1,:,:), objs(1,:), flags(1,:)] = gradient_optim(b.X_start_struct,p,b);	
x0s(1) = b.X_start_struct;

% 20 random ICs
parfor i = 2:num_runs
    [x0_vec,x0] = random_x0(b);
    [~, ~, ~, feasible_lin] = is_feasible(0, x0_vec, p, b);
    x0s(i) = x0;
    if feasible_lin
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
lin_feasible = isfinite(objs);
conv = flags > 0 & lin_feasible; % converged
kkt = flags == 1;
tol = 0.001;
opt = abs(objs - min(objs)) < tol; % this assumes that solutions converging 
% to the same objective value have the same optimal DVs, which is not guaranteed 
% but seems correct in this case by inspection of the table above

% optimal_and_converged = conv & opt;
% optimal_and_kkt = kkt & opt;
% 
% percent_lin_feas = sum(lin_feasible) / num_runs
% percent_converged_given_lin_feas = sum(conv) / sum(lin_feasible)
% percent_kkt_given_lin_feas = sum(kkt) / sum(lin_feasible)
% 
% percent_optimal_given_lin_feas = sum(opt) / num_runs
% percent_optimal_given_converged = sum(optimal_and_converged) ./ sum(conv)
% percent_optimal_given_kkt = sum(optimal_and_kkt) ./ sum(kkt)

edge_connections = [2 1;
                    3 1;
                    4 2; 
                    5 2;
                    6 2;
                    7 4;
                    8 4;
                    9 5; 
                    10 5;
                    11 6;
                    12 6];
node_names = {'1000 random points',...
    'linear feasible','linear infeasible',...
    'KKT','Converged, not KKT','Not converged',...
    'Global Min','Local Min','Global Min ','Local Min ','Global Min  ','Local Min  '};

edge_weights = [           sum(lin_feasible);                                                                               sum(~lin_feasible); ...
                sum(kkt);                      sum(conv & ~kkt);                      sum(~conv & lin_feasible);...
  sum(opt&kkt); sum(~opt&kkt);  sum(conv & ~kkt & opt); sum(conv & ~kkt & ~opt);  sum(~conv & opt); sum(~conv & ~opt & lin_feasible)];
                

G_J1 = digraph(edge_connections(:,2),edge_connections(:,1),edge_weights(:,1),node_names);
G_J2 = digraph(edge_connections(:,2),edge_connections(:,1),edge_weights(:,2),node_names);
figure
subplot(121)
plot(G_J1,'NodeLabel',G_J1.Nodes.Name,'EdgeLabel',G_J1.Edges.Weight,...
    'NodeFontSize',10,'EdgeFontSize',12,'MarkerSize',8,'LineWidth',1)
title('LCOE minimization')
subplot(122)
plot(G_J2,'NodeLabel',G_J2.Nodes.Name,'EdgeLabel',G_J2.Edges.Weight,...
    'NodeFontSize',10,'EdgeFontSize',12,'MarkerSize',8,'LineWidth',1)
title('Design Cost Minimization')
improvePlot

end
