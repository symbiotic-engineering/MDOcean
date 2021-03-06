function [Xs_opt, objs_opt, flags, probs] = gradient_optim(x0_input,p,b,which_objs)

if nargin == 0
    % set default parameters if function is run without input
    clc;close all
    p = parameters();
    b = var_bounds(p);
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
                                'MaxIterations',8);
                            
% iterate through material choices                            
for matl = 1%1:2:3 %b.M_min : b.M_max
    X = [D_f D_s_ratio h_f_ratio T_s_ratio F_max B_p w_n matl];

    [Xs_opt, objs_opt, flags, probs] = optimize_both_objectives(X,p,b,x0_input,opts,ploton,which_objs);

end

end

%%
function [Xs_opt, objs_opt, flags, probs] = optimize_both_objectives(X,p,b,x0_input,opts,ploton,which_objs)

    [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
       FOS_buckling, GM, P_elec, D_d, ~, g] = fcn2optimexpr(@simulation,X,p,...
                                            'OutputSize',{[1,1],[1,1],[1,2],...
                                            [1,1],[1,1],[1,1],[1,1],[1,1],....
                                            [1,1],[1,1],size(p.JPD),[1, 18]},...
                                            'ReuseEvaluation',true,'Analysis','off');%simulation(X, p);
    
    objs = [LCOE P_var];
    obj_names = {'LCOE','P_var'};
    probs = cell([1 2]);
    
    num_objectives = length(which_objs);
    Xs_opt = zeros(length(X),num_objectives);
    objs_opt = zeros(1,num_objectives);
    flags = zeros(1,num_objectives);

    % add constraints
    prob = optimproblem();
    prob.Constraints.Buoyancy_float_min     = B(1) >= 0;
    prob.Constraints.Buoyancy_float_max     = B(1) <= 1;
    prob.Constraints.Buoyancy_spar_min      = B(2) >= 0;
    prob.Constraints.Buoyancy_spar_max      = B(2) <= 1;
    prob.Constraints.FOS_float_hydro        = FOS1Y / p.FOS_min >= 1;
    prob.Constraints.FOS_column_hydro       = FOS2Y / p.FOS_min >= 1;
    prob.Constraints.FOS_plate_hydro        = FOS3Y / p.FOS_min >= 1;
    prob.Constraints.FOS_buckling_hydro     = FOS_buckling / p.FOS_min >= 1;
    prob.Constraints.GM                     = GM >= 0;
    prob.Constraints.P_positive             = P_elec >= 0;
    prob.Constraints.Damping                = D_d / p.D_d_min >= 1;
    prob.Constraints.Spar_height            = g(16) >= 0;
    prob.Constraints.LCOE_max               = g(17) >= 0;
    prob.Constraints.F_max_limit            = g(18) >= 0;

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
end
