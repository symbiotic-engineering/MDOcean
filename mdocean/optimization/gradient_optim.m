function [Xs_opt, objs_opt, flags, probs] = gradient_optim(x0_input,p,b,which_objs)

if nargin == 0
    % set default parameters if function is run without input
    clc;close all
    p = parameters();
    b = var_bounds();
    x0_input = b.X_start_struct;
    display = 'iter';
    plotfn = {@optimplotfval, @(x,~,~)optim_geomviz(x,p,b)};
    ploton = true;
else
    display = 'off';
    plotfn = [];
    ploton = false;
end

if nargin<4
    which_objs = [1 2]; % run both objectives by default
end

% create optimization variables for each of the design variables
sz = [1 1]; % create scalar variables
x1  = optimvar(b.var_names{1},  sz,'LowerBound',b.X_mins(1),  'UpperBound',b.X_maxs(1));
x2  = optimvar(b.var_names{2},  sz,'LowerBound',b.X_mins(2),  'UpperBound',b.X_maxs(2));
x3  = optimvar(b.var_names{3},  sz,'LowerBound',b.X_mins(3),  'UpperBound',b.X_maxs(3));
x4  = optimvar(b.var_names{4},  sz,'LowerBound',b.X_mins(4),  'UpperBound',b.X_maxs(4));
x5  = optimvar(b.var_names{5},  sz,'LowerBound',b.X_mins(5),  'UpperBound',b.X_maxs(5));
x6  = optimvar(b.var_names{6},  sz,'LowerBound',b.X_mins(6),  'UpperBound',b.X_maxs(6));
x7  = optimvar(b.var_names{7},  sz,'LowerBound',b.X_mins(7),  'UpperBound',b.X_maxs(7));
x8  = optimvar(b.var_names{8},  sz,'LowerBound',b.X_mins(8),  'UpperBound',b.X_maxs(8));
x9  = optimvar(b.var_names{9},  sz,'LowerBound',b.X_mins(9),  'UpperBound',b.X_maxs(9));
x10 = optimvar(b.var_names{10}, sz,'LowerBound',b.X_mins(10), 'UpperBound',b.X_maxs(10));
x11 = optimvar(b.var_names{11}, sz,'LowerBound',b.X_mins(11), 'UpperBound',b.X_maxs(11));
x12 = optimvar(b.var_names{12}, sz,'LowerBound',b.X_mins(12), 'UpperBound',b.X_maxs(12));
x13 = optimvar(b.var_names{13}, sz,'LowerBound',b.X_mins(13), 'UpperBound',b.X_maxs(13));
x14 = optimvar(b.var_names{14}, sz,'LowerBound',b.X_mins(14), 'UpperBound',b.X_maxs(14));

opts = optimoptions('fmincon',	'Display',display,...
                                'Algorithm','sqp',...%'interior-point',...
                                ...%"EnableFeasibilityMode",true,...
                                ...%"SubproblemAlgorithm","cg",...
                                'PlotFcn',plotfn,...
                                'MaxIterations',8,...
                                'FunValCheck','on',...
                                'ConstraintTolerance',1e-5);%,...
                                %'SpecifyObjectiveGradient',true,...
                                %'SpecifyConstraintGradient',true); % would require ALL constraints to be AD-supported

% iterate through material choices                            
for matl = 1%1:2:3 %b.M_min : b.M_max
    X = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 matl];
    [Xs_opt, objs_opt, flags, probs] = optimize_both_objectives(X,p,b,x0_input,opts,ploton,which_objs);

end

end

%%
function [Xs_opt, objs_opt, flags, probs] = optimize_both_objectives(X,p,b,x0_input,opts,ploton,which_objs)

    num_constraints = length(b.constraint_names);
    num_objectives = length(which_objs);

    [LCOE, P_var, ~, g] = fcn2optimexpr(@simulation,X,p,...
                                            'OutputSize',{[1,1],[1,1],size(p.JPD),[1, num_constraints]},...
                                            'ReuseEvaluation',true,'Analysis','off');%simulation(X, p);
    
    objs = [LCOE P_var];
    obj_names = {'LCOE','P_var'};
    probs = cell([1 2]); 
    
    Xs_opt = zeros(length(X),num_objectives);
    objs_opt = zeros(1,num_objectives);
    flags = zeros(1,num_objectives);

    % add constraints
    prob = optimproblem();
    for i = 1:num_constraints
        name = b.constraint_names{i};
        prob.Constraints.(name) = g(i) >= 0;
    end

    D_s   = X(1);     % inner diameter of float (m)
    D_f   = X(2);     % normalized diameter of spar column (-)
    T_f_2 = X(3);     % normalized draft of float (-)
    h_s   = X(4);     % normalized draft of spar (-)

    T_s = D_s * p.T_s_over_D_s;
    h_f = T_f_2 / p.T_f_2_over_h_f;
    D_d = p.D_d_over_D_s * D_s;

    MEEM = pi*p.harmonics / (p.besseli_argmax*2);
    prob.Constraints.linear_spar_natural_freq = D_d >= p.D_d_min;
    prob.Constraints.linear_float_spar_diam = D_s <= D_f - .01;
    prob.Constraints.linear_float_spar_draft = T_f_2 <= T_s - .01;
    prob.Constraints.linear_float_spar_tops = h_s - T_s >= h_f - T_f_2 + .01;
    prob.Constraints.linear_float_seafloor = p.h - T_f_2 >= MEEM * D_f; % M
    prob.Constraints.linear_spar_seafloor = p.h - T_s >= MEEM * D_s; % N

    % iterate through the two objectives: LCOE and P_var
    for i = 1:num_objectives
        which_obj = which_objs(i);
        prob.Objective = objs(which_obj);
        
        %show(prob)

        if length(x0_input)==1
            x0 = x0_input;
        elseif length(x0_input)==num_objectives
            x0 = x0_input(i);
        else
            error('x0 input struct has wrong size')
        end
            
        [X_opt_raw,obj_opt,flag,...
            output,lambda,grad,hess,problem] = run_solver(prob, obj_names{which_obj}, x0, opts, b.idxs_recover, b.filename_uuid);
        probs{i} = problem;

        tol = eps(2);
        min_active = abs(X_opt_raw - b.X_mins) < tol;
        max_active = abs(X_opt_raw - b.X_maxs) < tol;
        flexible_min_active = b.mins_flexible & min_active;
        flexible_max_active = b.maxs_flexible & max_active;
        if  any(flexible_min_active) || any(flexible_max_active)
            var_opt = b.var_names(1:end-1);
            max_string = strcat(var_opt(flexible_max_active)," ");
            min_string = strcat(var_opt(flexible_min_active)," ");
            msg = ['Optimization is up against a flexible variable bound, consider changing bounds. ',...
                   'At max: ' horzcat(max_string{:}) ', at min: ' horzcat(min_string{:})];
            warning(msg)
        end

        X_opt = [X_opt_raw; evaluate(X(end),struct())];   % add material back onto design vector
        [out(1),out(2)] = simulation(X_opt,p);          % rerun sim
        assert(out(which_obj) == obj_opt)               % check correct reordering of X_opt elements
        
        Xs_opt(:,i) = X_opt;
        objs_opt(i) = obj_opt;
        flags(i) = flag;

        % Post process
        if ploton
            plot_power_matrix(X_opt,p)
            visualize_geometry(X_opt,p)
        end
    end
    if ploton
        table_data = [Xs_opt(1:end-1,:), b.X_mins, b.X_maxs];
        objs_opt
        flags
        array2table(table_data,'RowNames',b.var_names(1:end-1),...
                'VariableNames',{'Min LCOE','Min cv','Min bound','Max bound'})
    end

end
