clear;clc;close all

p = parameters();
p.M = 1; % hardcode steel for now - run a for-loop
p.F_max = 1e6;

D_sft       = optimvar('D_sft',     [1 1],'LowerBound',0,   'UpperBound',50);
D_i_ratio   = optimvar('D_i_ratio', [1 1],'LowerBound',0,   'UpperBound',1);
D_or        = optimvar('D_or',      [1 1],'LowerBound',0,   'UpperBound',50);
N_WEC       = optimvar('N_WEC',     [1 1],'LowerBound',1,   'UpperBound',100); % integer but relaxed
D_int       = optimvar('D_int',     [1 1],'LowerBound',1e6, 'UpperBound',1e8);
mult        = optimvar('mult',   [length(p.Hs) length(p.T)],...
                                          'LowerBound',0,   'UpperBound',1);
slack1      = optimvar('slack1', [length(p.Hs) length(p.T)],'UpperBound',0);
slack2      = optimvar('slack2', [length(p.Hs) length(p.T)],'UpperBound',1);


X = [D_sft D_i_ratio D_or 0 N_WEC D_int];% K F_max];

[LCOE, D_env, Lt, B, FOS, P_elec, F_pt_unsat] = simulation(X, mult, p);

prob = optimproblem('Objective',LCOE);

opts = optimoptions('fmincon',	'Display','iter',...
                                        'Algorithm','sqp');

% these constraints, along with the upper bounds on slack1 and slack2 above,
% implement the constraint mult = min(p.F_max / F_pt, 1).
prob.Constraints.forceMultEquality      = mult + slack1 <= slack2;
prob.Constraints.forceMultActive        = slack2 <= p.F_max * F_pt_unsat.^-1;

show(prob)

%%
x0 = struct('D_sft',20,'D_i_ratio',.3,'D_or',30,'N_WEC',10,'D_int',1e7,...
    'mult',ones(size(p.JPD)),'slack1',zeros(size(p.JPD)),'slack2',ones(size(p.JPD)));
opt = solve(prob,x0,'Options',opts);
%%
X_opt = [opt.D_sft opt.D_i_ratio opt.D_or 0 opt.N_WEC opt.D_int];
[~,~,~,~,~,~,P_matrix,~] = simulation(X_opt, opt.mult, p);

[Hs,T] = meshgrid(p.T,p.Hs);
figure
contourf(Hs,T,P_matrix);
xlabel('Hs')
ylabel('T')
title('Power Matrix')
colorbar
