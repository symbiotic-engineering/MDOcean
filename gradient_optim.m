clear;clc;close all

p = parameters();
b = var_bounds(p);

D_sft       = optimvar('D_sft',     [1 1],'LowerBound',b.D_sft_min,    'UpperBound',b.D_sft_max);
D_i_ratio   = optimvar('D_i_ratio', [1 1],'LowerBound',b.D_i_ratio_min,'UpperBound',b.D_i_ratio_max);
D_or        = optimvar('D_or',      [1 1],'LowerBound',b.D_or_min,     'UpperBound',b.D_or_max);
N_WEC       = optimvar('N_WEC',     [1 1],'LowerBound',b.N_WEC_min,    'UpperBound',b.N_WEC_max);
D_int       = optimvar('D_int',     [1 1],'LowerBound',b.D_int_min,    'UpperBound',b.D_int_max);
w_n         = optimvar('w_n',       [1 1],'LowerBound',b.w_n_min,      'UpperBound',b.w_n_max);

x0 = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or',...
    b.D_or_nom,'N_WEC',b.N_WEC_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);

opts = optimoptions('fmincon',	'Display','iter',...
                                'Algorithm','sqp');
    
for matl = 1:2:3 %b.M_min : b.M_max
    X = [D_sft D_i_ratio D_or matl N_WEC D_int w_n];

    [LCOE, D_env, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, ~, ~] = fcn2optimexpr(@simulation,X,p);%simulation(X, p);

    prob = optimproblem('Objective',LCOE);
    
    % Structures Constraints
    prob.Constraints.B              = B >= p.B_min;
    prob.Constraints.FOS1YH         = FOS1Y(1) >= p.FOS_min;
    prob.Constraints.FOS1YP         = FOS1Y(2) >= p.FOS_min;
    prob.Constraints.FOS2YH         = FOS2Y(1) >= p.FOS_min;
    prob.Constraints.FOS2YP         = FOS2Y(2) >= p.FOS_min;
    prob.Constraints.FOS3YH         = FOS3Y(1) >= p.FOS_min;
    prob.Constraints.FOS3YP         = FOS3Y(2) >= p.FOS_min;
    prob.Constraints.FOS_bucklingH  = FOS_buckling(1) >= p.FOS_min;
    prob.Constraints.FOS_bucklingP  = FOS_buckling(2) >= p.FOS_min;
    prob.Constraints.GM             = GM >= 0;

    %show(prob)
    %% Run optimization
    [opt_x, opt_LCOE, ~,~,lambda] = solve(prob,x0,'Options',opts);
    
    %% Post process
    X_opt = [opt_x.D_sft opt_x.D_i_ratio opt_x.D_or matl opt_x.N_WEC opt_x.D_int opt_x.w_n];
    plot_power_matrix(X_opt,p)
end
