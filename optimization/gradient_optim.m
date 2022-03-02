function [Xs_opt, objs_opt, flag, problem] = gradient_optim(x0_input,p,b)

if nargin == 0
    % set default parameters if function is run without input
    clc;close all
    p = parameters();
    b = var_bounds(p);
    x0_input = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or_ratio',...
        b.D_or_ratio_nom,'N_WEC',b.N_WEC_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);
    display = 'iter';
    plotfn = @optimplotfval;
    ploton = true;
else
    display = 'off';
    plotfn = [];
    ploton = false;
end

% create optimization variables for each of the design variables
D_sft       = optimvar('D_sft',     [1 1],'LowerBound',b.D_sft_min,     'UpperBound',b.D_sft_max);
D_i_ratio   = optimvar('D_i_ratio', [1 1],'LowerBound',b.D_i_ratio_min, 'UpperBound',b.D_i_ratio_max);
D_or_ratio  = optimvar('D_or_ratio',[1 1],'LowerBound',b.D_or_ratio_min,'UpperBound',b.D_or_ratio_max);
N_WEC       = optimvar('N_WEC',     [1 1],'LowerBound',b.N_WEC_min,     'UpperBound',b.N_WEC_max);
D_int       = optimvar('D_int',     [1 1],'LowerBound',b.D_int_min,     'UpperBound',b.D_int_max);
w_n         = optimvar('w_n',       [1 1],'LowerBound',b.w_n_min,       'UpperBound',b.w_n_max);

opts = optimoptions('fmincon',	'Display',display,...
                                'Algorithm','sqp',...
                                'PlotFcn',plotfn,...
                                'MaxIterations',10);
                            
% iterate through material choices                            
for matl = 1%1:2:3 %b.M_min : b.M_max
    X = [D_sft D_i_ratio D_or_ratio matl N_WEC D_int w_n];

    [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, ~] = fcn2optimexpr(@simulation,X,p);%simulation(X, p);
    
    prob1 = optimproblem('Objective',LCOE);
    prob2 = optimproblem('Objective',P_var);
    objs = {'LCOE','P_var'};
    probs = {prob1 prob2};
    
    num_objectives = 2;
    Xs_opt = zeros(length(X),num_objectives);
    objs_opt = zeros(1,num_objectives);

    % iterate through the two objectives: LCOE and P_var
    for i=1:num_objectives
        prob = probs{i};
        % add constraints
        prob.Constraints.Buoyancy               = B/p.B_min >= 1;
        prob.Constraints.FOS_float_hydro        = FOS1Y(1)/p.FOS_min >= 1;
        prob.Constraints.FOS_float_ptrain       = FOS1Y(2)/p.FOS_min >= 1;
        prob.Constraints.FOS_column_hydro       = FOS2Y(1)/p.FOS_min >= 1;
        prob.Constraints.FOS_column_ptrain      = FOS2Y(2)/p.FOS_min >= 1;
        prob.Constraints.FOS_plate_hydro        = FOS3Y(1)/p.FOS_min >= 1;
        prob.Constraints.FOS_plate_ptrain       = FOS3Y(2)/p.FOS_min >= 1;
        prob.Constraints.FOS_buckling_hydro     = FOS_buckling(1)/p.FOS_min >= 1;
        prob.Constraints.FOS_buckling_ptrain    = FOS_buckling(2)/p.FOS_min >= 1;
        prob.Constraints.GM                     = GM >= 0;
        prob.Constraints.P_positive             = P_elec >= 0;

        %show(prob)

        if length(x0_input)==1
            x0 = x0_input;
        elseif length(x0_input)==num_objectives
            x0 = x0_input(i);
        else
            error('x0 input struct has wrong size')
        end
        
        solver_based = true;
        %% Run optimization
        % create folder for generated objectives if it doesn't already exist        
        if solver_based
            generated_folder = 'optimization/generated';
            if ~exist(generated_folder,'dir')
                mkdir(generated_folder)
                addpath(generated_folder)
            end
            % Convert to solver-based
            problem = prob2struct(prob,x0,...
                'ObjectiveFunctionName',['generatedObjective' objs{i}],...
                'FileLocation',generated_folder);
            problem.options = opts;
            
            % Run unscaled optimization
            [X_opt,obj_opt,flag,~,~,~,hess_unscaled] = fmincon(problem);
            
            if flag <= 0 % if the unscaled optimization did not arrive at an optimal
                % use the unscaled optimization hessian to find scale factor
                scale = 1./sqrt(diag(hess_unscaled));

                % Formulate a new scaled optimization problem     
                problem_s = problem;
                problem_s.options.MaxIterations = 100;
                problem_s.objective = @(x) problem.objective(x .* scale);  
                problem_s.nonlcon   = @(x) problem.nonlcon(x .* scale);

                inv_scale = 1./(scale);
                problem_s.x0 = inv_scale .* X_opt;
                problem_s.lb = inv_scale .* problem.lb;
                problem_s.ub = inv_scale .* problem.ub;

                % Run scaled optimization problem
                [X_opt,obj_opt,flag,output,lambda,grad,hess] = fmincon(problem_s);
                X_opt = scale .* X_opt;
            end
            
            % Rearrange and check outputs
            X_opt = [X_opt(4); X_opt(1); X_opt(3); matl; X_opt(5); X_opt(2); X_opt(6)]; % reorder elements based on order in autogenerated objective files
            [out(1),out(2)] = simulation(X_opt,p);
            assert(out(i) == obj_opt) % check that reordering of X_opt elements is correct
            
        else
            [opt_x, obj_opt, flag,output,lambda] = solve(prob,x0,'Options',opts);
            X_opt = [opt_x.D_sft opt_x.D_i_ratio opt_x.D_or_ratio matl opt_x.N_WEC opt_x.D_int opt_x.w_n];
        end
        
        Xs_opt(:,i) = X_opt;
        objs_opt(i) = obj_opt;

        %% Post process
        if ploton
            plot_power_matrix(X_opt,p)
            visualize_geometry(X_opt,p)
        end
    end

end
