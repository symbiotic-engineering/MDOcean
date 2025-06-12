function [x,fval,timeEpsConstraint,timeGradientFree,output] = pareto_search(filename_uuid,num_seeds_vector)
    if nargin<1
        filename_uuid='';
    end
    if nargin<2
        num_seeds_vector = 20;
    end
    
    p = parameters();
    b = var_bounds();
    b.filename_uuid = filename_uuid;
    num_DVs = length(b.X_starts);
    
    turnLCOEtoPower = true;
    objFcnName = 'generatedObjective';
    objFcn1 = str2func([objFcnName b.obj_names{1}]);
    objFcn2 = str2func([objFcnName b.obj_names{2}]);

    %% Calculate seed points for the pareto front
    
    t = tic;
    [X0,fvals,probs] = get_seeds_epsilon_constraint(p,b,num_DVs,objFcn2,turnLCOEtoPower,max(num_seeds_vector));
    timeEpsConstraint = toc(t)
    disp('Finished finding pareto seed points. Now starting paretosearch.')

    %% Gradient free search to fill in pareto front
    t2 = tic;
    [x,fval,flag,output,residuals,probMO] = run_gradient_free(p,X0,fvals,probs,num_DVs,objFcn1,objFcn2,turnLCOEtoPower,solver);
    timeGradientFree = toc(t2)

    percentSeed = timeEpsConstraint / (timeEpsConstraint + timeGradientFree) * 100

    %% Process and save results
    utopia = min(fval);
    hold on
    
    plot(utopia(1),utopia(2),'gp','MarkerFaceColor','g','MarkerSize',20)

    % show which constraints are active along the pareto front
    tol = probMO.options.ConstraintTolerance;
    idx = constraint_active_plot(residuals,fval,tol,b);

    x_sorted = x(idx,b.idxs_recover)

    % save mat file to be read by pareto_heuristics.m
    date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
    save(['optimization/multiobjective/pareto_search_results_' date '.mat'],"fval","x","residuals","tol","p")
end

function [probMO, X0_struct] = create_multiobj_prob(X0,fvals,scale,probs,num_DVs,objFcn1,objFcn2,turnLCOEtoPower)
    probMO = probs{1};
    
    if turnLCOEtoPower
        p_zero_design_cost = p;
        % hack that makes cost constant, so minimizing LCOE is actually maximizing average power
        p_zero_design_cost.cost_perN_mult = 0;
        p_zero_design_cost.cost_perW_mult = 0;
        p_zero_design_cost.cost_perkg_mult = [0 0 0];

        objFcn1Revised = @(x)objFcn1new(x,objFcn1,p_zero_design_cost);
        scale(1) = 1;
    else
        objFcn1Revised = @(x)objFcn1(x,{p});
    end

    probMO.objective = @(x)[objFcn1Revised(x)*scale(1), objFcn2(x,{p})*scale(2)];
    
    X0_struct = struct('X0',X0,'Fvals',fvals .* scale);
    probMO.nvars = num_DVs;
    if turnLCOEtoPower
        X0_struct.Fvals = [];
    end
end

function [x,fval,flag,output,residuals,probMO] = run_gradient_free(p,X0,fvals,probs,num_DVs,objFcn1,objFcn2,turnLCOEtoPower,solver)
    scale = [10 1];
    [probMO_base, X0_struct] = create_multiobj_prob(X0,fvals,scale,probs,num_DVs,objFcn1,objFcn2,turnLCOEtoPower);

    if strcmpi(solver,'patternsearch')
        [x,fval,flag,output,residuals,probMO] = pattern_search(probMO_base, X0_struct);
    elseif strcmpi(solver,'genetic')
        error('GA not yet implemented')
    else
        error('invalid solver')
    end

    if flag==-2
        error(output.message)
    end

    if turnLCOEtoPower
        for i=1:length(fval)
            fval(i,1) = objFcn1(x(i,:),{p}); % put LCOE back in to avoid confusing pareto curve heuristics
        end
    end

    fval = fval ./ repmat(scale,length(fval),1);
end

function [x,fval,flag,output,residuals,probMO] = pattern_search(probMO, X0_struct)
    %% Set up pareto search algorithm

    probMO.options = optimoptions('paretosearch','Display','iter',...
        'PlotFcn','psplotparetof','MinPollFraction',1,...
        'ParetoSetChangeTolerance',1.6e-8,'MaxIterations',100,...
        'UseParallel',true);
    if ~isempty(X0_struct.X0)
        probMO.options.InitialPoints = X0_struct;
    else
        probMO.x0 = [];
        warning('All starting points are infeasible')
    end
    probMO.solver = 'paretosearch';

    %% Execute pareto search
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

end

function [X0,fvals,probs] = get_seeds_epsilon_constraint(p,b,num_DVs,objFcn2,turnLCOEtoPower,num_seeds_vector)
    x0 = b.X_start_struct;
    idxs = b.idxs_sort;

    % get pareto front endpoints by running optimization
    [X_opt,objs_opt,flag_opt,probs] = gradient_optim(x0,p,b);
    [~,P_var_min_LCOE] = simulation(X_opt(:,1),p);
    
    LCOE_max = p.LCOE_max;
    LCOE_min = objs_opt(1);
    P_var_min = objs_opt(2);
    fvals_opt = [LCOE_min P_var_min_LCOE; LCOE_max P_var_min];

    % max power seed
    if turnLCOEtoPower
        [X_opt_pwr,val,flag_max_pwr] = max_avg_power(p,b);
        X_opt(:,end+1) = X_opt_pwr;
        fvals_opt(end+1,:) = [val.LCOE val.J_capex_design];
        flag_opt(end+1) = flag_max_pwr;
    end

    % remove infeasible single-objective optima
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

    [X_seeds,LCOE_seeds,P_var_seeds,init_failed] = get_seeds_scalar(max(num_seeds_vector));

    % remove failed seeds
    X_seeds(init_failed,:) = [];
    P_var_seeds(init_failed) = [];
    LCOE_seeds(init_failed) = [];
    

    % combine seeds
    X0 = [X_0; X_opt(idxs,:)'; X_seeds];
    fvals = [fvals_0; fvals_opt; LCOE_seeds' P_var_seeds'];
end

function [X_seeds,LCOE_seeds,P_var_seeds] = get_seeds_scalar(num_seeds,LCOE_min,LCOE_max,num_DVs,p,b)
    theta = linspace(0,pi/2,num_seeds+2);
    theta = theta(2:end-1);
    LCOE_seeds = LCOE_min + (1-sin(theta))*(LCOE_max-LCOE_min);
    X_seeds = zeros(length(LCOE_seeds),num_DVs);
    P_var_seeds = zeros(1,length(LCOE_seeds));
    init_failed = false(1,length(LCOE_seeds));

    parfor i = 1:length(LCOE_seeds)
        new_p = p;
        new_p.LCOE_max = LCOE_seeds(i);
        which_obj = 2;
        new_b = b;
        if isempty(new_b.filename_uuid)
            t = getCurrentTask();
            if ~isempty(t)
                new_b.filename_uuid = num2str(t.ID);
            end
        end
        [X_opt_tmp,obj_tmp,flag_tmp] = gradient_optim(x0,new_p,new_b,which_obj);
        if flag_tmp == -2 % Initial pareto point is not feasible - this prevents a LPalg error in paretosearch
            init_failed(i) = true;
            warning('Initial pareto point not feasible (fmincon returned -2 flag), removing.')
        else
            X_seeds(i,:) = X_opt_tmp(idxs)';
    
            % debugging checks on optimization convergence and objective values
            obj_check = objFcn2(X_opt_tmp(idxs)',{new_p});
            assert(obj_tmp == obj_check)
            [~, P_var_seeds(i)] = simulation(X_opt_tmp, new_p);
            assert(obj_tmp == P_var_seeds(i))
        end
    end
end

function neg_pwr_per_cost = objFcn1new(x,oldFcnLCOE,p_zero_design_cost)

   % cost0 refers to the part of the numerator of LCOE that does not depend 
   % on design: cost0 = capex0 + opex0
   cost0_per_power = oldFcnLCOE(x,{p_zero_design_cost});
   pwr_per_cost0 = 1/cost0_per_power;
   neg_pwr_per_cost = -pwr_per_cost0;

   clear generatedFunction_simulation1_withReuse % required since using different parameters for the two objs
   
end
