function [x,fval,lambda] = pareto_search(filename_uuid)
    if nargin==0; filename_uuid=''; end
    
    p = parameters();
    b = var_bounds();
    
    b.filename_uuid = filename_uuid;
    x0 = b.X_start_struct;
    idxs = b.idxs_sort;
    num_DVs = length(b.X_starts);
    
    turnLCOEtoPower = true;
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
    num_seeds = 20;
    theta = linspace(0,pi/2,num_seeds+2);
    theta = theta(2:end-1);
    LCOE_seeds = LCOE_min + (1-sin(theta))*(LCOE_max-LCOE_min);
    X_seeds = zeros(length(LCOE_seeds),num_DVs);
    P_var_seeds = zeros(1,length(LCOE_seeds));
    init_failed = false(1,length(LCOE_seeds));

    for i = 1:length(LCOE_seeds)
        new_p = p;
        new_p.LCOE_max = LCOE_seeds(i);
        which_obj = 2;
        new_b = b;
        if isempty(new_b.filename_uuid)
            t = getCurrentTask();
            if ~isempty(t)
                %new_b.filename_uuid = num2str(t.ID);
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

    % remove failed seeds
    X_seeds(init_failed,:) = [];
    P_var_seeds(init_failed) = [];
    LCOE_seeds(init_failed) = [];

    %% Set up pareto search algorithm
    probMO = probs{1};
    scale = [10 1];
    
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
    
    X0 = [X_0; X_opt(idxs,:)'; X_seeds];
    fvals = [fvals_0; fvals_opt; LCOE_seeds' P_var_seeds'];
    X0_struct = struct('X0',X0,'Fvals',fvals .* scale);

    if turnLCOEtoPower
        X0_struct.Fvals = [];
    end

    probMO.options = optimoptions('paretosearch','Display','iter',...
        'PlotFcn','psplotparetof','MinPollFraction',1,...
        'ParetoSetChangeTolerance',1.6e-8,'MaxIterations',100,...
        'UseParallel',true);
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

    if turnLCOEtoPower
        for i=1:length(fval)
            fval(i,1) = objFcn1(x(i,:),{p}); % put LCOE back in to avoid confusing pareto curve heuristics
        end
    end

    fval = fval ./ repmat(scale,length(fval),1);

    % show which constaints are active along the pareto front
    tol = probMO.options.ConstraintTolerance;
    idx = constraint_active_plot(residuals,fval,tol);

    x_sorted = x(idx,b.idxs_recover);
% solve for lambda values
    [rows,~] = size(x);
    for i = 1:length(LCOE_seeds)
        p.LCOE_max = LCOE_seeds(i);
        which_objs = 2;
    end

    for i = 1:rows
        input = x_sorted(i,:);
        %x_0_new.(b.var_names{1}) = input(1);
        %x_0_new.(b.var_names{2}) = input(2);
        %x_0_new.(b.var_names{3}) = input(3);
        %x_0_new.(b.var_names{4}) = input(4);
        %x_0_new.(b.var_names{5}) = input(5);
        %x_0_new.(b.var_names{6}) = input(6);
        %x_0_new.(b.var_names{7}) = input(7);
        %x_0_new.(b.var_names{8}) = input(8);
        %x_0_new.(b.var_names{9}) = input(9);
        %x_0_new.(b.var_names{10}) = input(10);
        %x_0_new.(b.var_names{11}) = input(11);
        %x_0_new.(b.var_names{12}) = input(12);
        %x_0_new.(b.var_names{13}) = input(13);
        %x_0_new.(b.var_names{14}) = input(14);

        x_0_new = cell2struct(num2cell(input)',b.var_names(1:end-1),1);

        [Xs_opt_original, ~, ~, ~, lambda_original, gs] = gradient_optim(x_0_new,p,b,which_objs);
        g(:,i) = gs;
        Xs_opt_original(8,:) = [];
        Xs_opt(:,i) = Xs_opt_original;
        lambda.active(:,i) = lambda_original.ineqnonlin;
        lambda.lower(:,i) = lambda_original.lower;
        lambda.upper(:,i) = lambda_original.upper;
    end
    residuals2.ineqnonlin = g';
    residuals2.lower = Xs_opt'-b.X_mins';
    residuals2.upper = Xs_opt'-b.X_maxs';
    [~]=constraint_active_plot(residuals2,fval,tol);

    % save mat file to be read by pareto_heuristics.m
    date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
    save(['optimization/multiobjective/pareto_search_results_' date '.mat'],"fval","x","residuals","tol","p","lambda")
end

function neg_pwr_per_cost = objFcn1new(x,oldFcnLCOE,p_zero_design_cost)

   % cost0 refers to the part of the numerator of LCOE that does not depend 
   % on design: cost0 = capex0 + opex0
   cost0_per_power = oldFcnLCOE(x,{p_zero_design_cost});
   pwr_per_cost0 = 1/cost0_per_power;
   neg_pwr_per_cost = -pwr_per_cost0;

   clear generatedFunction_simulation1_withReuse % required since using different parameters for the two objs
   
end

function [idx] = constraint_active_plot(residuals,fval,tol)
    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    con_active = abs(residuals.ineqnonlin) < tol;
    [~,idx] = sort(fval(:,1)); % order by increasing LCOE

    figure
    subplot 311

    %subplot 311
    H=subplot(3,1,1);
    spy(lb_active(idx,:)');
    title('Lower Bound Active')
    yticks(linspace(1,8,8));
    yticklabels({'WEC Surface Float Outer Diameter', ...
        'Ratio of WEC Surface Float Inner Diameter to Outer Diameter', ...
        'Ratio of WEC Surface Float Height to Outer Diameter', ...
        'Percent of WEC Spar Submergence','Maximum Powertrain Force', ...
        'Powertrain/Controller Damping','Controller Natural Frequency', ...
        'Material'});
    xlabel("Designs, ordered by increasing LCOE", "Fontsize", 9)
    grid on

    subplot 312
    %subplot 312
    I=subplot(3,1,2);
    spy(ub_active(idx,:)')
    title('Upper Bound Active')
    yticks(linspace(1,8,8))
    yticklabels({'WEC Surface Float Outer Diameter', ...
        'Ratio of WEC Surface Float Inner Diameter to Outer Diameter', ...
        'Ratio of WEC Surface Float Height to Outer Diameter', ...
        'Percent of WEC Spar Submergence','Maximum Powertrain Force', ...
        'Powertrain/Controller Damping','Controller Natural Frequency', ...
        'Material'});
    xlabel("Designs, ordered by increasing LCOE", "Fontsize", 9)
    grid on

    subplot 313
    J=subplot(3,1,3);
    spy(con_active(idx,:)')
    title('Constraint Active')
    xlabel("Designs, ordered by increasing LCOE", "Fontsize", 9)
    yticks(linspace(1,14,14));
    yticklabels({'Prevent Float Too Heavy', 'Prevent Float Too Light', ...
        'Prevent Spar Too Heavy','Prevent Spar Too Light', ...
        'Metacentric Height','Float Yield','Column Yield','Reaction Plate Yield' ...
        'Column Buckling','Net Generated Power','Minimum Damping Plate Diameter', ...
        'Prevent Float Above Top of the Spar','Prevent LCOE Greater Than Nominal', ...
        'Prevent Irrelevant Max Force'});
    grid on

    improvePlot

    %subplot 311
    set(H, "FontSize",10)
    title(H,'Lower Bound Active', "FontSize", 16)

    %subplot 312
    set(I, "FontSize",10)
    title(I,'Upper Bound Active', "FontSize", 16)

    %subplot 313
    set(J, "FontSize",9.5)
    title(J,'Constraint Active', "FontSize", 16)

end