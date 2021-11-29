function [X_opt, opt_LCOE, flag] = gradient_optim(x0,p,b)

if nargin == 0
    clc;close all
    p = parameters();
    b = var_bounds(p);
    x0 = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or',...
        b.D_or_nom,'N_WEC',b.N_WEC_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);
    display = 'iter';
    plotfn = @optimplotfval;
    ploton = true;
else
    display = 'off';
    plotfn = [];
    ploton = false;
end

D_sft       = optimvar('D_sft',     [1 1],'LowerBound',b.D_sft_min,    'UpperBound',b.D_sft_max);
D_i_ratio   = optimvar('D_i_ratio', [1 1],'LowerBound',b.D_i_ratio_min,'UpperBound',b.D_i_ratio_max);
D_or        = optimvar('D_or',      [1 1],'LowerBound',b.D_or_min,     'UpperBound',b.D_or_max);
N_WEC       = optimvar('N_WEC',     [1 1],'LowerBound',b.N_WEC_min,    'UpperBound',b.N_WEC_max);
D_int       = optimvar('D_int',     [1 1],'LowerBound',b.D_int_min,    'UpperBound',b.D_int_max);
w_n         = optimvar('w_n',       [1 1],'LowerBound',b.w_n_min,      'UpperBound',b.w_n_max);

opts = optimoptions('fmincon',	'Display',display,...
                                'Algorithm','sqp',...
                                'PlotFcn',plotfn);
                            
for matl = 1%1:2:3 %b.M_min : b.M_max
    X = [D_sft D_i_ratio D_or matl N_WEC D_int w_n];

    [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, ~] = fcn2optimexpr(@simulation,X,p);%simulation(X, p);

    prob = optimproblem('Objective',LCOE);

    % Structures Constraints
    prob.Constraints.Buoyancy               = B/p.B_min >= 1;
    prob.Constraints.FOS_float_hydro        = FOS1Y(1)/p.FOS_min >= 1;
    prob.Constraints.FOS_float_ptrain       = FOS1Y(2)/p.FOS_min >= 1;
    prob.Constraints.FOS_column_hydro       = FOS2Y(1)/p.FOS_min >= 1;
    prob.Constraints.FOS_column_ptrain      = FOS2Y(2)/p.FOS_min >= 1;
    prob.Constraints.FOS_plate_hydro        = FOS3Y(1)/p.FOS_min >= 1;
    prob.Constraints.FOS_plate_ptrain       = FOS3Y(2)/p.FOS_min >= 1;
    prob.Constraints.FOS_buckling_hydro     = FOS_buckling(1)/p.FOS_min >= 1;
    prob.Constraints.FOS_buckling_ptrain    = FOS_buckling(2)/p.FOS_min >= 1;
    prob.Constraints.GM             = GM >= 0;
    prob.Constraints.P_positive     = P_elec >= 0;

    %show(prob)
    
    solver_based = true;
    %% Run optimization
    if solver_based
        % Convert to solver-based
        problem = prob2struct(prob,x0);
        problem.options = opts;
        [X_opt_un,opt_LCOE,flag,output,lambda,grad,hess_un] = fmincon(problem)
                Diag=diag(hess_un)
                cond(hess_un)
                %% Scaled Optimization
                Scale=1./sqrt(Diag); %Scale Factor Calculation
                problem_s=problem
                problem_s.objective=@(x)generatedObjective(x.*Scale, {p});
                problem_s.nonlcon=@(x)generatedConstraints(x.*Scale, {p,p,p,p,p,p,p,p,p,p,p});
                X_sx0=1./(Scale);
                problem_s.x0=X_sx0.*X_opt_un;
                problem_s.lb=X_sx0.*problem.lb;
                problem_s.ub=X_sx0.*problem.ub;
                %X_opts=X_scale.*X_opt %Scaled design variables
                [X_opt,opt_LCOE,flag,output,lambda,grad,hess] = fmincon(problem_s)
                Diag_s=diag(hess)
             X_opt = [X_opt(1) X_opt(2) X_opt(3) matl X_opt(4) X_opt(5) X_opt(6)]   
             
    else
    	[opt_x, opt_LCOE, flag,output,lambda] = solve(prob,x0,'Options',opts);
        X_opt = [opt_x.D_sft opt_x.D_i_ratio opt_x.D_or matl opt_x.N_WEC opt_x.D_int opt_x.w_n];
    end
    
    %% Post process
    if ploton
        plot_power_matrix(Scale.*X_opt,p)
    end
end

end
