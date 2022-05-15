clc;clear
p = parameters();
b = var_bounds(p);
x0 = b.X_start_struct;

[X_opt,objs_min,~, probs] = gradient_optim(x0,p,b);

%%
probMO = probs{1};
probMO.objective = @(x)[generatedObjectiveLCOE(x,{p}), generatedObjectiveP_var(x,{p})];

idxs = [1 6 2 5 4 3 7];
X0 = [probMO.x0'; X_opt(idxs,:)'];
probMO.options = optimoptions('paretosearch','Display','iter','PlotFcn','psplotparetof',...
    'InitialPoints',X0,'MinPollFraction',.9);%'ParetoSetSize',100);
probMO.solver = 'paretosearch';
probMO.nvars = 7;
%%
[x,fval] = paretosearch(probMO)
utopia = min(fval);
hold on

plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)
