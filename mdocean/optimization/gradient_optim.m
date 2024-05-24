function [Xs_opt, objs_opt, flags, probs] = gradient_optim(x0_input,p,b,which_objs,materials)

if nargin == 0
    % set default parameters if function is run without input
    clc;close all
    p = parameters();
    b = var_bounds();
    x0_input = b.X_start_struct;
    display = 'iter';
    plotfn = @optimplotfval;
    ploton = true;
else
    display = 'off';
    plotfn = [];
    ploton = false;
end

if nargin<4
    which_objs = [1 2]; % run both objectives by default
end

if nargin<5
    materials = [1 3]; % run both steels by default
    % materials = b.M_min : b.M_max
end

% create optimization variables for each of the design variables
sz = [1 1]; % create scalar variables
D_f         = optimvar('D_f',       sz,'LowerBound',b.D_f_min,       'UpperBound',b.D_f_max);
D_s_ratio   = optimvar('D_s_ratio', sz,'LowerBound',b.D_s_ratio_min, 'UpperBound',b.D_s_ratio_max);
h_f_ratio   = optimvar('h_f_ratio', sz,'LowerBound',b.h_f_ratio_min, 'UpperBound',b.h_f_ratio_max);
T_s_ratio   = optimvar('T_s_ratio', sz,'LowerBound',b.T_s_ratio_min, 'UpperBound',b.T_s_ratio_max);
F_max       = optimvar('F_max',     sz,'LowerBound',b.F_max_min,     'UpperBound',b.F_max_max);
B_p         = optimvar('D_int',     sz,'LowerBound',b.B_p_min,       'UpperBound',b.B_p_max);
w_n         = optimvar('w_n',       sz,'LowerBound',b.w_n_min,       'UpperBound',b.w_n_max);

opts = optimoptions('fmincon',	'Display',display,...
                                'Algorithm','sqp',...
                                'PlotFcn',plotfn,...
                                'MaxIterations',8,...
                                'FunValCheck','on',...
                                'ConstraintTolerance',1e-5);

num_matls = length(materials);
num_objs = length(which_objs);
num_DVs = length(b.var_names);

% initialize loop outputs
Xs_opt_tmp = zeros(num_DVs,num_objs,num_matls);
objs_opt_tmp = zeros(num_objs,num_matls);
flags_tmp = zeros(num_objs,num_matls);
probs_tmp = cell(num_objs,num_matls);

% iterate through material choices
for matl_idx = 1 : num_matls
    matl = materials(matl_idx);
    X = [D_f D_s_ratio h_f_ratio T_s_ratio F_max B_p w_n matl];

    [Xs_opt_tmp(:,:,matl_idx), ...
     objs_opt_tmp(:,matl_idx), ...
     flags_tmp(:,matl_idx), ...
     probs_tmp(:,matl_idx)] = optimize_both_objectives(X,p,b,x0_input,opts,ploton,which_objs);

end

% post process to find the optimal material for each objective
opt_matl_idx = zeros(1, num_objs);

% if the optimization only succeeded for one material, use that material
success = flags_tmp >= 0;
num_successful_matls_each_obj = sum(success, 2);
one_success = num_successful_matls_each_obj == 1;
no_success = num_successful_matls_each_obj == 0;
multiple_success = ~one_success & ~no_success;
[~,opt_matl_idx(one_success)] = find(success(one_success,:));
% if the optimization did not succeed for any material, use the material with the lowest objectives
[~,opt_matl_idx(no_success)] = min(objs_opt_tmp(no_success,:));
% if the optimization succeeded for multiple materials, choose from among the
% successful materials and select the material with the lower objective
multiple_success_objs = objs_opt_tmp(multiple_success,:);
multiple_success_objs(~success(multiple_success,:)) = NaN;
[~,opt_matl_idx(multiple_success)] = min(multiple_success_objs,[],2);

% collect outputs at the optimal material - requires strange indexing
row_idx = repmat((1:num_DVs)', [1 num_objs]);
col_idx = repmat(1:num_objs, [num_DVs 1]);
page_idx = repmat(opt_matl_idx, [num_DVs 1]);
X_opt_idx = sub2ind(size(Xs_opt_tmp), row_idx, col_idx, page_idx);
Xs_opt   = Xs_opt_tmp(X_opt_idx);
two_dim_idx = sub2ind(size(objs_opt_tmp), 1:num_objs, opt_matl_idx);
objs_opt = objs_opt_tmp(two_dim_idx);
flags    = flags_tmp(two_dim_idx);
probs    = probs_tmp(two_dim_idx);

end

%%
function [Xs_opt, objs_opt, flags, probs] = optimize_both_objectives(X,p,b,x0_input,opts,ploton,which_objs)

    num_constraints = length(b.constraint_names);
    num_objectives = length(which_objs);

    [LCOE, P_var, ~, g] = fcn2optimexpr(@simulation,X,p,...
                                            'OutputSize',{[1,1],[1,1],size(p.JPD),[1, num_constraints]},...
                                            'ReuseEvaluation',true,'Analysis','off');%simulation(X, p);
    
    objs = [LCOE P_var];
    obj_names = {'LCOE','P_var'};
    probs = cell([1 num_objectives]); 
    
    Xs_opt = zeros(length(X),num_objectives);
    objs_opt = zeros(1,num_objectives);
    flags = zeros(1,num_objectives);

    % add constraints
    prob = optimproblem();
    for i = 1:num_constraints
        name = b.constraint_names{i};
        prob.Constraints.(name) = g(i) >= 0;
    end

    % iterate through the two objectives: LCOE and P_var
    for i = 1:num_objectives
        which_obj = which_objs(i);
        prob.Objective = objs(which_obj);
        
        %show(prob)

        if length(x0_input)==1
            x0 = x0_input;
        elseif length(x0_input)==num_objectives
            x0 = x0_input(i);
        else
            error('x0 input struct has wrong size')
        end
        
        [X_opt_raw,obj_opt,flag,output,lambda,grad,hess,problem] = run_solver(prob, obj_names{which_obj}, x0, opts);
        probs{i} = problem;

                       % D_f   D_s_ratio h_f_ratio T_s_ratio F_max B_p w_n]
        mins_flexible = [false false     false     false     false true  true]';
        maxs_flexible = [true  false     false     false     true  true  true]';
        tol = eps(2);
        if any(abs(X_opt_raw(mins_flexible) - b.X_mins(mins_flexible)) < tol) ...
                || any(abs(X_opt_raw(maxs_flexible) - b.X_maxs(maxs_flexible)) < tol)
            warning('Optimization is up against a flexible variable bound, consider changing bounds')
        end

        X_opt = [X_opt_raw; evaluate(X(8),struct())];   % add material back onto design vector
        [out(1),out(2)] = simulation(X_opt,p);          % rerun sim
        assert(out(which_obj) == obj_opt)               % check correct reordering of X_opt elements
        
        Xs_opt(:,i) = X_opt;
        objs_opt(i) = obj_opt;
        flags(i) = flag;

        % Post process
        if ploton
            plot_power_matrix(X_opt,p)
            visualize_geometry(X_opt,p)
        end
    end
    if ploton
        table_data = [Xs_opt(1:end-1,:), b.X_mins, b.X_maxs];
        objs_opt
        flags
        array2table(table_data,'RowNames',b.var_names(1:end-1),...
                'VariableNames',{'Min LCOE','Min cv','Min bound','Max bound'})
    end

end
