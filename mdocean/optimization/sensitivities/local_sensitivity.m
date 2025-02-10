
b = var_bounds();
[p,T] = parameters();

%param_names = T.name_pretty(T.sweep);  % list of parameters to sweep
%params = T.name(T.sweep);
params = {'h'};

dvar_names = b.var_names_pretty;

%%
% use the optimal x as x0 to speed up the sweeps
% and obtain gradients
x0 = b.X_start_struct;
[x0_vec, J0, ~, ~, lambdas, grads, hesses] = gradient_optim(x0,p,b);

num_constr = length(b.constraint_names);

disp('done optimizing, starting param sensitivity')
tic
local_sens_both_obj_all_param(x0_vec, p, params, lambdas, grads, hesses, num_constr)
toc

function [LCOE_L,  X_LCOE_L, ...
          P_var_L, X_Pvar_L] = local_sens_deltas(x0, p, param_name, ratios, matrix, lambdas, J0)

    delta_p = p0 * (ratios - 1);

    [par_x_star_par_p, dJdp, dJstar_dp] = get_local_sens();

    % deltas based on specific ratios
    delta_x_star = par_x_star_par_p * delta_p;
    delta_J = dJdp * delta_p; 
    delta_Jstar = dJstar_dp * delta_p;

    dxs(:,obj) = delta_x_star;
    dJs(obj)   = delta_Jstar;

    x_new = x0s + dxs;
    J_new = J0 + dJs;
end

function [par_x_star_par_p, dJstar_dp, ...
            dJdp, par_lam_par_p] = local_sens_both_obj_all_param(x0s, p, params, lambdas, grads, hesses, num_constr)
    for obj = 1:2
        x0 = x0s(:,obj);
        lambda = lambdas(obj);
        grad = grads(:,obj);
        hess = hesses(:,:,obj);

        [par_x_star_par_p, dJstar_dp, ...
            dJdp, par_lam_par_p] = local_sens_one_obj_all_param(x0, p, params, lambda, grad, hess, num_constr, obj);
    end
end

function [par_x_star_par_p, dJstar_dp, ...
            dJdp, par_lam_par_p] = local_sens_one_obj_all_param(x0, p, params, lambda, grad, hess, num_constr, obj)
% construct matrix equation for a single J - slide 47 of lec 10 of MDO

    % lagrange multpliers
    lambda_nl  = lambda.ineqnonlin;
    lambda_lin = lambda.ineqlin;
    lambda_lb = lambda.lower;
    lambda_ub = lambda.upper;

    % active constraints
    active     = lambda_nl ~= 0;
    active_lin = lambda_lin ~= 0;
    active_lb  = lambda_lb ~= 0;
    active_ub  = lambda_ub ~= 0;
    % fixme: need to actually use these bounds

    % left hand side matrix - independent of parameter
    matrix = local_sens_LHS_matrix(x0,p,hess,active,active_lin,active_lb,active_ub);
    
    % sweep parameters
    for i = 1:length(params)
        param_name = params{i};
        % right hand side vector - depends on parameter
        [vector, par_J_par_p, dJdp] = local_sens_RHS_vector_and_dJdp(obj, p, param_name, x0, num_constr, ...
                                                        lambda_nl, lambda_lin, lambda_lb, lambda_ub, ...
                                                        active, active_lin, active_lb, active_ub);
    
        % matrix solution
        sol = matrix \ vector;
        par_x_star_par_p = sol(1:length(x0)-1);
        par_lam_par_p = sol(length(x0):end);
    
        % total derivative
        par_J_par_x = grad;
        dJstar_dp = par_J_par_p + par_J_par_x.' * par_x_star_par_p;

        %delta_p_nl_becomes_inactive = -lambda_nl(active) ./ par_lam_par_p
        %delta_p_becomes_active
    end
end

function [vector, par_J_par_p, dJdp] = local_sens_RHS_vector_and_dJdp(obj, p, param_name, x0, num_constr, ...
                                                        lambda_nl, lambda_lin, lambda_lb, lambda_ub, ...
                                                        active, active_lin, active_lb, active_ub)


    p0 = p.(param_name);
    [par_J_par_p, par_g_par_p, ...
        par_par_J_par_x_par_p, ...
        par_par_g_par_x_par_p,...
        par_g_lin_par_p,...
        par_par_g_lin_par_x_par_p] = get_partials(obj, p, param_name, x0, p0,num_constr);

    % assume that x bounds are parameter-independent
    par_g_lb_par_p = zeros(size(active_lb)); par_g_ub_par_p = zeros(size(active_ub));

    % use all constraints for J sensitivity - MDO book - just J, not J*
    dJdp = par_J_par_p ...
    - lambda_nl.' * par_g_par_p ...
    - lambda_lin.'* par_g_lin_par_p ...
    - lambda_lb.' * par_g_lb_par_p ...
    - lambda_ub.' * par_g_ub_par_p;

    % only use active constraints for J* and x* sensitivity
    % d vector
    d = [par_g_par_p(active);
        par_g_lin_par_p(active_lin);
        par_g_lb_par_p(active_lb);
        par_g_ub_par_p(active_ub)];

    % c vector
    c_obj_term = par_par_J_par_x_par_p;
    c_nl_constraint_term = par_par_g_par_x_par_p * lambda_nl;
    c_lin_constraint_term = par_par_g_lin_par_x_par_p * lambda_lin;
    c_lb_term = 0; %par_par_g_lb_par_x_par_p * lambda_lb;
    c_ub_term = 0; %par_par_g_ub_par_x_par_p * lambda_ub;
    c = c_obj_term + c_nl_constraint_term + c_lin_constraint_term + c_lb_term + c_ub_term;
    vector = -[c; d];
end

function matrix = local_sens_LHS_matrix(x0,p,hess,active,active_lin,active_lb,active_ub)

    % A is just the Hessian of the Lagrangian, which we already have from the optimization
    A = hess;

    % number of active constraints
    num_constr_active = sum(active);
    num_constr_active_tot = sum([active; active_lin; active_lb; active_ub]);

    % obtain B = dg/dx with finite difference for nonlinear constraints
    active_constraint_handle = @(X) get_constraints([X;1],p,active);
    B_nl = finite_difference_vector_x(active_constraint_handle,x0(1:end-1),num_constr_active);

    % obtain B from constraint A_ineq matrix for linear constraints and bounds
    [~, A_ineq_lin_part] = is_feasible(0, x0, p);
    A_ineq_lin = zeros(size(A_ineq_lin_part,1), length(x0)-1);
    A_ineq_lin(:, 1:size(A_ineq_lin_part,2)) = A_ineq_lin_part;
    B_lin = A_ineq_lin(active_lin,:).';
    A_ineq_lb = diag(-active_lb);
    A_ineq_ub = diag(active_ub);
    B_lb = -A_ineq_lb(:,active_lb);
    B_ub = A_ineq_ub(:,active_ub);
    B = [B_nl B_lin B_lb B_ub];

    % make matrix
    matrix = [A B; B.' zeros(num_constr_active_tot)];
end

function [par_J_par_p, par_g_par_p, ...
            par_par_J_par_x_par_p, ...
            par_par_g_par_x_par_p,...
            par_g_lin_par_p,...
            par_par_g_lin_par_x_par_p] = get_partials(obj, p, param_name, x0, p0,num_constr)

        %  compute par_y_par_p using finite difference
        y_fcn_handle = @(param_value) get_sim_outputs_and_derivs(obj,p,param_name,param_value,x0);
        par_y_par_p = finite_difference_scalar_x(y_fcn_handle,p0);
    
        % decompose par_y_par_p into its sections
        par_J_par_p = par_y_par_p(1,1);
        par_g_par_p = par_y_par_p(1,2 : 1+num_constr).';
        par_par_J_par_x_par_p = par_y_par_p(2:end, 1);
        par_par_g_par_x_par_p = par_y_par_p(2:end, 2 : 1+num_constr);

        % linear sensitivities (analytical)
        [~, ~, dAdp, dbdp] = lin_ineq_constraints(p, param_name);
        num_constr_lin = size(dAdp,1);
        num_DVs_constr_lin = size(dAdp,2);
        par_g_lin_par_p = dAdp * x0(1:num_DVs_constr_lin) - dbdp; % zeros(1,num_constr_lin);
        par_par_g_lin_par_x_par_p = zeros(length(x0)-1, num_constr_lin);
        par_par_g_lin_par_x_par_p(1:num_DVs_constr_lin,:) = dAdp.';
end

function y = get_sim_outputs_and_derivs(obj,p,param_name,param_value,x0)
    p.(param_name) = param_value;

    % F = [J(obj),g]
    F_handle = @(X)get_sim_outputs([X;1],p,obj);

    % get F directly from simulation
    F = F_handle(x0);

    % get parF_parX from finite difference
    parF_parx = finite_difference_vector_x(F_handle,x0(1:end-1),length(F));

    % assemble y from F and parF_parX
    y = [F;parF_parx];
end

% wrappers to get the right outputs to make function handles
function F = get_sim_outputs(X,p,obj)
    [J1, J2, ~, g] = simulation(X, p);
    J = [J1,J2];
    F = [J(obj),-g]; % g negative because sim is g>0, sens wants g<0
end
function g_active = get_constraints(X,p,active)
    [~, ~, ~, g] = simulation(X, p);
    g_active = -g(active); % g negative because sim is g>0, sens wants g<0
end

function deriv = finite_difference_vector_x(fcn_handle,x_1_vec,ny)
% fcn handle: a function expecting input vector x and returning output matrix y (size ny)
    
    nx = length(x_1_vec);
    deriv = zeros([nx ny]);
    for i=1:nx
        x_before = x_1_vec(1:i-1);
        x_after  = x_1_vec(i+1:end);
        scalar_handle = @(x_i) fcn_handle([x_before; x_i; x_after]);
        deriv(i,:) = finite_difference_scalar_x(scalar_handle,x_1_vec(i));
    end
end

function deriv = finite_difference_scalar_x(f_handle,x_1)
% fcn handle: a function expecting input scalar x and returning output matrix y
    assert(isscalar(x_1))
    y_1 = f_handle(x_1);
    delta_x = max(1e-10, x_1*2e-4);
    x_2 = x_1 + delta_x;
    y_2 = f_handle(x_2);

    delta_y = y_2 - y_1;
    deriv = delta_y / delta_x;
    % deriv is size(y)
end

