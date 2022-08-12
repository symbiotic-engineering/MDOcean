function [x,fval] = pareto_search()
    p = parameters();
    b = var_bounds(p);
    x0 = b.X_start_struct;
    
    %% Calculate seed points for the pareto front
    % get pareto front endpoints by running optimization
    [X_opt,objs_opt,~,probs] = gradient_optim(x0,p,b);

    % get extra seed values in the middle by optimizing at fixed LCOEs
    LCOE_max = p.LCOE_max;
    LCOE_seeds = [0.1 0.2 0.3 0.4];
    X_seeds = zeros(length(LCOE_seeds),7);
    P_var_seeds = zeros(1,length(LCOE_seeds));
    for i = 1:length(LCOE_seeds)
        p.LCOE_max = LCOE_seeds(i);
        which_obj = 2;
        [X_opt_tmp,obj_tmp,flag_tmp] = gradient_optim(x0,p,b,which_obj);
        idxs = [1 6 2 5 4 3 7];
        X_seeds(i,:) = X_opt_tmp(idxs)';

        % debugging checks on optimization convergence and objective values
        assert(flag_tmp==1);
        obj_check = generatedObjectiveP_var(X_opt_tmp(idxs)',{p});
        assert(obj_tmp == obj_check)
        [~, P_var_seeds(i)] = simulation(X_opt_tmp, p);
        assert(obj_tmp == P_var_seeds(i))
    end
    p.LCOE_max = LCOE_max;

    %% Set up pareto search algorithm
    probMO = probs{1};
    scale = [10 0.1];
    probMO.objective = @(x)[generatedObjectiveLCOE(x,{p})*scale(1), generatedObjectiveP_var(x,{p})*scale(2)];
    
    [LCOE_x0,P_var_x0] = simulation([b.X_starts; 1],p);
    [~,P_var_min_LCOE] = simulation(X_opt(:,1),p);
    X0 = [probMO.x0'; X_opt(idxs,:)'; X_seeds];
    fvals = [LCOE_x0,P_var_x0; objs_opt(1) P_var_min_LCOE; LCOE_max objs_opt(2); LCOE_seeds' P_var_seeds'];
    X0_struct = struct('X0',X0,'Fvals',fvals .* scale);

    probMO.options = optimoptions('paretosearch','Display','iter',...
        'PlotFcn','psplotparetof','InitialPoints',X0_struct,'MinPollFraction',1,...
        'ParetoSetChangeTolerance',1.6e-8,'MaxIterations',100);
    probMO.solver = 'paretosearch';
    probMO.nvars = 7;
    %% Execute pareto search
    [x,fval,flag,output,residuals] = paretosearch(probMO);

    %% Process and save results
    utopia = min(fval);
    hold on
    
    plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)

    fval = fval ./ repmat(scale,length(fval),1);

    % show which constaints are active along the pareto front
    tol = probMO.options.ConstraintTolerance;
    constraint_active_plot(residuals,fval,tol)

    % save mat file to be read by pareto_bruteforce.m
    save('optimization/multiobjective/pareto_search_results6',"fval","x","residuals")
end

function [] = constraint_active_plot(residuals,fval,tol)
    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    con_active = abs(residuals.ineqnonlin) < tol;

    [~,idx] = sort(fval(:,1)); % order by increasing LCOE

    figure
    subplot 311
    spy(lb_active(idx,:)');
    title('Lower Bound Active')

    subplot 312
    spy(ub_active(idx,:)')
    title('Upper Bound Active')

    subplot 313
    spy(con_active(idx,:)')
    title('Constraint Active')
end