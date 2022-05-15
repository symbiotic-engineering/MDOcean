clc
p = parameters();
b = var_bounds(p);
x0 = b.X_start_struct;

[x_opt,objs_min,~, prob] = gradient_optim(x0,p,b);
idxs = 1:7;
X_opt = x_opt(idxs);
%%
prob.Objective.LCOE = @(x)generatedObjectiveLCOE(x,{p});
prob.Objective.Pvar = @(x)generatedObjectiveP_var(x,{p});
prob.options = optimoptions('paretosearch','Display','iter','PlotFcn','psplotparetof',...
    'InitialPoints',[prob.x0'; X_opt'],'MinPollFraction',.9);%'ParetoSetSize',100);
prob.solver = 'paretosearch';
prob.nvars = 7;
%%
[x,fval] = paretosearch(prob);
utopia = min(fval);
hold on

plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)
