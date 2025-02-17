function [x,fval] = pareto_search(filename_uuid)
    if nargin==0; filename_uuid=''; end
    
    p = parameters();
    b = var_bounds();
    
    b.filename_uuid = filename_uuid;
    x0 = b.X_start_struct;
    idxs = b.idxs_sort;
    num_DVs = length(b.X_starts);
    
    objFcnName = 'generatedObjective';
    objFcn1 = str2func([objFcnName b.obj_names{1}]);
    objFcn2 = str2func([objFcnName b.obj_names{2}]);

    %% Calculate seed points for the pareto front
    % get pareto front endpoints by running optimization
    [X_opt,objs_opt,flag_opt,probs] = gradient_optim(x0,p,b);
    [~,P_var_min_LCOE] = simulation(X_opt(:,1),p);
    
    LCOE_max = p.LCOE_max;
    LCOE_min = objs_opt(1);
    P_var_min = objs_opt(2);

    % remove infeasible single-objective optima
    fvals_opt = [LCOE_min P_var_min_LCOE; LCOE_max P_var_min];
    single_obj_failed = flag_opt == -2;
    X_opt(:,single_obj_failed) = [];
    fvals_opt(single_obj_failed,:) = [];

    % use x0 start as seed, if x0 is feasible
    [LCOE_x0,P_var_x0,~,g_0] = simulation([b.X_starts; 1],p);
    x0_feasible = is_feasible(g_0, [b.X_starts; 1], p, b);
    if x0_feasible
        fvals_0 = [LCOE_x0,P_var_x0];
        X_0 = b.X_starts(idxs);
    else
        fvals_0 = [];
        X_0 = [];
    end

    % get extra seed values in the middle by optimizing at fixed LCOEs
    if LCOE_max < LCOE_min
        LCOE_max = 2*LCOE_min;
    end
    num_seeds = 8;
    LCOE_seeds = linspace(LCOE_min, LCOE_max, num_seeds+2);
    LCOE_seeds = LCOE_seeds(2:end-1); % remove min and max, since we already have those from gradient optim
    LCOE_seeds = [LCOE_seeds(1) (LCOE_seeds(1)+LCOE_seeds(2))/2 LCOE_seeds(2:end)]; % add extra seed between first and second to fill gap there
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
            obj_check = objFcn2(X_opt_tmp(idxs)',{p});
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
    
    probMO.objective = @(x)[objFcn1(x,{p})*scale(1), objFcn2(x,{p})*scale(2)];
    
    X0 = [X_0; X_opt(idxs,:)'; X_seeds];
    fvals = [fvals_0; fvals_opt; LCOE_seeds' P_var_seeds'];
    X0_struct = struct('X0',X0,'Fvals',fvals .* scale);

    probMO.options = optimoptions('paretosearch','Display','iter',...
        'PlotFcn','psplotparetof','MinPollFraction',1,...
        'ParetoSetChangeTolerance',1.6e-8,'MaxIterations',100);
    if ~isempty(X0)
        probMO.options.InitialPoints = X0_struct;
    else
        probMO.x0 = [];
        warning('All starting points are infeasible')
    end
    probMO.solver = 'paretosearch';
    probMO.nvars = num_DVs;
    %% Execute pareto search
    disp('Finished finding pareto seed points. Now starting paretosearch.')
    try
        [x,fval,flag,output,residuals] = paretosearch(probMO);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:emptyObjectDotAssignment'))
            msg = ['Pareto search could not find any feasible starting points, ' ...
                   'so plot errored. Retrying without plot.'];
            warning(msg)
            probMO.options.PlotFcn = [];
            [x,fval,flag,output,residuals] = paretosearch(probMO);
        else
            rethrow(ME)
        end
    end
    if flag==-2
        error(output.message)
    end

    %% Process and save results
    utopia = min(fval);
    hold on
    
    plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)

    fval = fval ./ repmat(scale,length(fval),1);

    % show which constaints are active along the pareto front
    tol = probMO.options.ConstraintTolerance;
    idx = constraint_active_plot(residuals,fval,tol,b);

    x_sorted = x(idx,b.idxs_recover)

    % save mat file to be read by pareto_heuristics.m
    date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
    save(['optimization/multiobjective/pareto_search_results_' date '.mat'],"fval","x","residuals","tol")
end
