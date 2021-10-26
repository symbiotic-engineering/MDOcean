                                
D_sft       = optimvar('D_sft',     [1 1],'LowerBound',0,   'UpperBound',50);
D_i_ratio   = optimvar('D_i_ratio', [1 1],'LowerBound',0,   'UpperBound',1);
D_or        = optimvar('D_or',      [1 1],'LowerBound',0,   'UpperBound',50);
N_WEC       = optimvar('N_WEC',     [1 1],'LowerBound',1,   'UpperBound',100, 'Type','integer');
D_int       = optimvar('D_int',     [1 1],'LowerBound',1e6, 'UpperBound',1e8);

X = [D_sft D_i_ratio D_or material N_WEC D_int];% K F_max];

p = parameters();
p.M = 1; % hardcode steel for now - run a for-loop

prob = optimproblem('Objective',simulation(X, p));
prob.x0 = [20 .3 30 10 1e7];
prob.options = optimoptions('fmincon',	'Display','iter',...
                                        'Algorithm','sqp');

prob.Constraints.c1 = D_sft > 1; % fixme with real constraints

show(prob)

%%
X = solve(prob,x0);

