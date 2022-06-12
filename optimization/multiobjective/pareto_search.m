function [x,fval] = pareto_search()
    p = parameters();
    b = var_bounds(p);
    x0 = b.X_start_struct;
    
    [X_opt,~,~, probs] = gradient_optim(x0,p,b);
    %%
    LCOE = p.LCOE_max;
    LCOE_seeds = [0.1 0.2 0.35 0.5 0.7];
    X_seeds = zeros(length(LCOE_seeds),7);
    for i = 1:length(LCOE_seeds)
        p.LCOE_max = LCOE_seeds(i);
        [X_opt_tmp,~,~,~] = gradient_optim(x0,p,b);
        idxs = [1 6 2 5 4 3 7];
        X_seeds(i,:) = X_opt_tmp(idxs,2)';
    end
    p.LCOE_max = LCOE;
    
    %%
    probMO = probs{1};
    scale = [10 0.1];
    probMO.objective = @(x)[generatedObjectiveLCOE(x,{p})*scale(1), generatedObjectiveP_var(x,{p})*scale(2)];
    
    
    X0 = [probMO.x0'; X_opt(idxs,:)'; X_seeds];
    probMO.options = optimoptions('paretosearch','Display','iter','PlotFcn','psplotparetof',...
        'InitialPoints',X0, 'MinPollFraction',1,'ParetoSetChangeTolerance',1.6e-8,'MaxIterations',100);
    probMO.solver = 'paretosearch';
    probMO.nvars = 7;
    %%
    [x,fval,flag,output,residuals] = paretosearch(probMO);
    utopia = min(fval);
    hold on
    
    plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)

    fval = fval ./ repmat(scale,length(fval),1);

    tol = probMO.options.ConstraintTolerance;
    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    con_active = abs(residuals.ineqnonlin) < tol;

    [~,idx] = sort(fval(:,1)); % order by increasing LCOE

    figure
    subplot 311
    spy(lb_active(idx,:)');
    subplot 312
    spy(ub_active(idx,:)')
    subplot 313
    spy(con_active(idx,:)')


    save('optimization/multiobjective/pareto_search_results6',"fval","x","residuals")
end