clc
p = parameters();
b = var_bounds(p);
x0 = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or',...
        b.D_or_nom,'N_WEC',b.N_WEC_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);

[x_opt,LCOE_min,~, prob] = gradient_optim(x0,p,b);
%idxswap = [2 6 3 1 5 7];
idxswap = [1:3 5:7];
X_opt = x_opt(idxswap);

prob.objective = @(x)[generatedObjectiveLCOE(x,{p}) generatedObjectiveP_var(x,{p})];
prob.options = optimoptions('paretosearch','Display','iter','PlotFcn','psplotparetof',...
    'InitialPoints',[prob.x0'; X_opt'],'MinPollFraction',.9);%'ParetoSetSize',100);
prob.solver = 'paretosearch';
prob.nvars = 6;
%%
[x,fval] = paretosearch(prob);
utopia = min(fval);
hold on

plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)
