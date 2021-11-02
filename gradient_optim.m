clear;clc;close all

p = parameters();
matl = 1; % hardcode steel for now - run a for-loop

w_max = 2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
w_min = 2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
D_sft       = optimvar('D_sft',     [1 1],'LowerBound',0,    'UpperBound',50);
D_i_ratio   = optimvar('D_i_ratio', [1 1],'LowerBound',0,    'UpperBound',1);
D_or        = optimvar('D_or',      [1 1],'LowerBound',0,    'UpperBound',50);
N_WEC       = optimvar('N_WEC',     [1 1],'LowerBound',1,    'UpperBound',100); % integer but relaxed
D_int       = optimvar('D_int',     [1 1],'LowerBound',1e6,  'UpperBound',1e8);
w_n         = optimvar('w_n',       [1 1],'LowerBound',w_min,'UpperBound',w_max);

X = [D_sft D_i_ratio D_or matl N_WEC D_int w_n];

[LCOE, D_env, B, FOS1Y, FOS2Y, ...
            FOS3Y, FOS_buckling,GM...
            ~, ~, F_pt_unsat] = fcn2optimexpr(@simulation, X, p);

prob = optimproblem('Objective',LCOE);

opts = optimoptions('fmincon',	'Display','iter',...
                                        'Algorithm','sqp');
 %Structures Constraints
prob.Constraints.B              = B >= p.B_min;
prob.Constraints.FOS1YH         = FOS1Y(1) >= p.FOS_min;
prob.Constraints.FOS1YP         = FOS1Y(2) >= p.FOS_min;
prob.Constraints.FOS2YH         = FOS2Y(1) >= p.FOS_min;
prob.Constraints.FOS2YP         = FOS2Y(2) >= p.FOS_min;
prob.Constraints.FOS3YH         = FOS3Y(1) >= p.FOS_min;
prob.Constraints.FOS3YP         = FOS3Y(2) >= p.FOS_min;
prob.Constraints.FOS_bucklingH  = FOS_buckling(1) >= p.FOS_min;
prob.Constraints.FOS_bucklingP  = FOS_buckling(2) >= p.FOS_min; 
prob.Constraints.GM             = GM <=0;
% these constraints, along with the upper bounds on slack1 and slack2 above,
% implement the constraint mult = min(p.F_max / F_pt, 1).
% prob.Constraints.forceMultEquality      = mult + slack1 <= slack2;
% prob.Constraints.forceMultActive        = slack2 <= p.F_max * F_pt_unsat.^-1;

show(prob)

%%
x0 = struct('D_sft',20,'D_i_ratio',.3,'D_or',30,'N_WEC',10,'D_int',1e7,'w_n',2*pi/8);
    %'mult',ones(size(p.JPD)),'slack1',zeros(size(p.JPD)),'slack2',ones(size(p.JPD)));
opt = solve(prob,x0,'Options',opts);
%%
X_opt = [opt.D_sft opt.D_i_ratio opt.D_or matl opt.N_WEC opt.D_int opt.w_n];

plot_power_matrix(X_opt,p)
