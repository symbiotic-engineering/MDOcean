function [x,fval] = pareto_search(filename_uuid)
    if nargin==0; filename_uuid=''; end
    
    p = parameters();
    b = var_bounds();
    b.filename_uuid = filename_uuid;
    x0 = b.X_start_struct;
    idxs = b.idxs_sort;
    num_DVs = length(b.X_starts);
    
    %% Calculate seed points for the pareto front
    % get pareto front endpoints by running optimization
    [X_opt,objs_opt,~,probs] = gradient_optim(x0,p,b);

    % get extra seed values in the middle by optimizing at fixed LCOEs
    LCOE_max = p.LCOE_max;
    LCOE_min = objs_opt(1);
    P_var_min = objs_opt(2);

    num_seeds = 8;
    LCOE_seeds = linspace(LCOE_min, LCOE_max, num_seeds+2);
    LCOE_seeds = LCOE_seeds(2:end-1); % remove min and max, since we already have those from gradient optim
    X_seeds = zeros(length(LCOE_seeds),num_DVs);
    P_var_seeds = zeros(1,length(LCOE_seeds));
    init_failed = false(1,length(LCOE_seeds));
    
    for i = 1:length(LCOE_seeds)
        p.LCOE_max = LCOE_seeds(i);
        which_obj = 2;
        [X_opt_tmp,obj_tmp,flag_tmp] = gradient_optim(x0,p,b,which_obj);
        if flag_tmp == -2 % Initial pareto point is not feasible - this prevents a LPalg error in paretosearch
            init_failed(i) = true;
            warning('Initial pareto point not feasible (fmincon returned -2 flag), removing.')
        else
            X_seeds(i,:) = X_opt_tmp(idxs)';
    
            % debugging checks on optimization convergence and objective values
            obj_check = generatedObjectiveP_var(X_opt_tmp(idxs)',{p});
            assert(obj_tmp == obj_check)
            [~, P_var_seeds(i)] = simulation(X_opt_tmp, p);
            assert(obj_tmp == P_var_seeds(i))
        end
    end

    % remove failed seeds
    X_seeds(init_failed,:) = [];
    P_var_seeds(init_failed) = [];
    LCOE_seeds(init_failed) = [];

    p.LCOE_max = LCOE_max;

    %% Set up pareto search algorithm
    probMO = probs{1};
    scale = [10 0.1];
    probMO.objective = @(x)[generatedObjectiveLCOE(x,{p})*scale(1), generatedObjectiveP_var(x,{p})*scale(2)];
    
    [LCOE_x0,P_var_x0] = simulation([b.X_starts; 1],p);
    [~,P_var_min_LCOE] = simulation(X_opt(:,1),p);
    X0 = [probMO.x0'; X_opt(idxs,:)'; X_seeds];
    fvals = [LCOE_x0,P_var_x0; LCOE_min P_var_min_LCOE; LCOE_max P_var_min; LCOE_seeds' P_var_seeds'];
    X0_struct = struct('X0',X0,'Fvals',fvals .* scale);

    probMO.options = optimoptions('paretosearch','Display','iter',...
        'PlotFcn','psplotparetof','InitialPoints',X0_struct,'MinPollFraction',1,...
        'ParetoSetChangeTolerance',1.6e-8,'MaxIterations',100);
    probMO.solver = 'paretosearch';
    probMO.nvars = num_DVs;
    %% Execute pareto search
    disp('Finished finding pareto seed points. Now starting paretosearch.')
    [x,fval,flag,output,residuals] = paretosearch(probMO);

    %% Process and save results
    utopia = min(fval);
    hold on
    
    plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)

    fval = fval ./ repmat(scale,length(fval),1);

    % show which constaints are active along the pareto front
    tol = probMO.options.ConstraintTolerance;
    idx = constraint_active_plot(residuals,fval,tol);

    x_sorted = x(idx,b.idxs_recover);

    % save mat file to be read by pareto_heuristics.m
    date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
    save(['optimization/multiobjective/pareto_search_results_' date '.mat'],"fval","x","residuals")
end

function [idx] = constraint_active_plot(residuals,fval,tol)
    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    nlcon_active = abs(residuals.ineqnonlin) < tol;
    lincon_active = abs(residuals.ineqlin) < tol;

    % merge sea state slamming constraints
    nlcon_active(:,16) = any(nlcon_active(:,16:end),2);
    nlcon_active(:,17:end) = [];

    [~,idx] = sort(fval(:,1)); % order by increasing LCOE

    figure
    subplot 221
    spy(lb_active(idx,:)');
    title('Lower Bound Active')

    subplot 222
    spy(ub_active(idx,:)')
    title('Upper Bound Active')

    subplot 223
    spy(nlcon_active(idx,:)')
    title('Nonlinear Constraint Active')

    subplot 224
    spy(lincon_active(idx,:)')
    title('Linear Constraint Active')
end