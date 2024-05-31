function [x,fval] = pareto_search()
    p = parameters();
    b = var_bounds();
    x0 = b.X_start_struct;
    
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
    X_seeds = zeros(length(LCOE_seeds),7);
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
            idxs = [1 6 2 5 4 3 7];
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
    idx = constraint_active_plot(residuals,fval,tol);

    cols = [1 3 6 5 4 2 7];
    x_sorted = x(idx,cols);

    % solve for lambda values
    [rows,~] = size(x);
    for i = 1:length(LCOE_seeds)
        p.LCOE_max = LCOE_seeds(i);
        which_objs = 2;
    end

    for i = 1:rows
        input = x(i,:);
        x_0_new.D_f = input(1);
        x_0_new.D_s_ratio = input(2);
        x_0_new.h_f_ratio = input(3);
        x_0_new.T_s_ratio = input(4);
        x_0_new.F_max = input(5);
        x_0_new.D_int = input(6);
        x_0_new.w_n = input(7);
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

    % save mat file to be read by pareto_bruteforce.m
    save('optimization/multiobjective/pareto_search_results',"fval","x","residuals",...
    "lambda")

end

function [idx] = constraint_active_plot(residuals,fval,tol)
    lb_active = abs(residuals.lower) < tol;
    ub_active = abs(residuals.upper) < tol;
    con_active = abs(residuals.ineqnonlin) < tol;

    [~,idx] = sort(fval(:,1)); % order by increasing LCOE

    figure

    %subplot 311
    H=subplot(3,1,1, 'Position', [0.195 0.787 0.86 0.16]);
    spy(lb_active(idx,:)');
    yticks(linspace(1,8,8));
    yticklabels({'WEC Surface Float Outer Diameter', ...
        'Ratio of WEC Surface Float Inner Diameter to Outer Diameter', ...
        'Ratio of WEC Surface Float Height to Outer Diameter', ...
        'Percent of WEC Spar Submergence','Maximum Powertrain Force', ...
        'Powertrain/Controller Damping','Controller Natural Frequency', ...
        'Material'});
    xlabel("Designs, ordered by increasing LCOE", "Fontsize", 9)
    grid on

    %subplot 312
    I=subplot(3,1,2, 'Position', [0.195 0.505 0.86 0.16]);
    spy(ub_active(idx,:)')
    yticks(linspace(1,8,8))
    yticklabels({'WEC Surface Float Outer Diameter', ...
        'Ratio of WEC Surface Float Inner Diameter to Outer Diameter', ...
        'Ratio of WEC Surface Float Height to Outer Diameter', ...
        'Percent of WEC Spar Submergence','Maximum Powertrain Force', ...
        'Powertrain/Controller Damping','Controller Natural Frequency', ...
        'Material'});
    xlabel("Designs, ordered by increasing LCOE", "Fontsize", 9)
    grid on

    %subplot 313
    J=subplot(3,1,3, 'Position', [0.2 0.08 0.85 0.3]);
    spy(con_active(idx,:)')
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