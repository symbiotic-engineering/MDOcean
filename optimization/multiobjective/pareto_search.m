function [x,fval] = pareto_search()
    p = parameters();
    b = var_bounds(p);
    x0 = b.X_start_struct;
    
    [X_opt,~,~, probs] = gradient_optim(x0,p,b);
    %%
    LCOE = p.LCOE_max;
    p.LCOE_max = 0.2;
    [X_opt_2,~,~,~] = gradient_optim(x0,p,b);
    p.LCOE_max = 0.1;
    [X_opt_3,~,~,~] = gradient_optim(x0,p,b);
    p.LCOE_max = LCOE;
    
    %%
    probMO = probs{1};
    probMO.objective = @(x)[generatedObjectiveLCOE(x,{p}), generatedObjectiveP_var(x,{p})];
    
    idxs = [1 6 2 5 4 3 7];
    X0 = [probMO.x0'; X_opt(idxs,:)'; X_opt_2(idxs,2)'; X_opt_3(idxs,2)'];
    probMO.options = optimoptions('paretosearch','Display','iter','PlotFcn','psplotparetof',...
        'InitialPoints',X0, 'MinPollFraction',1,'ParetoSetChangeTolerance',2e-8);
    probMO.solver = 'paretosearch';
    probMO.nvars = 7;
    %%
    [x,fval,flag,output] = paretosearch(probMO);
    utopia = min(fval);
    hold on
    
    plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)
end