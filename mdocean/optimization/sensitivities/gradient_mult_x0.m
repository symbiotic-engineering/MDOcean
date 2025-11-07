function [X_opt,objs,flags,x0s] = gradient_mult_x0(p,b)

num_runs = 1000;
num_DVs = length(b.var_names);
num_objs = 2;
objs = Inf(num_runs,num_objs);
X_opt = NaN(num_runs,num_DVs,num_objs);
flags = zeros(num_runs,num_objs);

% nominal ICs
[X_opt(1,:,:), objs(1,:), flags(1,:)] = gradient_optim(b.X_start_struct,p,b);
x0s(1) = b.X_start_struct;

% many random ICs
parfor i = 2:num_runs
    [x0_vec,x0] = random_x0(b);
    [~, ~, ~, feasible_lin] = is_feasible(0, x0_vec, p, b);
    x0s(i) = x0;
    if feasible_lin
        [X_opt(i,:,:), objs(i,:), flags(i,:)] = gradient_optim(x0,p,b);
    end
end

end
