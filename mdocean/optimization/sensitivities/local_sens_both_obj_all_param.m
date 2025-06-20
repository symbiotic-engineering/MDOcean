function [par_x_star_par_p_norm, dJstar_dp_norm, ...
            dJdp_norm, par_J_par_p_norm, delta_p_norm] = local_sens_both_obj_all_param(x0s, J0, p, params, p_val, param_idxs, lambdas, grads, hesses, num_constr)
    for obj = 1%:2
        x0 = x0s(:,obj);
        lambda = lambdas(obj);
        grad = grads(:,obj);
        hess = hesses(:,:,obj);

        [par_x_star_par_p, dJstar_dp, ...
            dJdp, par_J_par_p, ...
            delta_p_change_activity] = local_sens_one_obj_all_param(x0, p, params, param_idxs, lambda, grad, hess, num_constr, obj);
        
        % normalization
        dJstar_dp_norm   = dJstar_dp.'   .* p_val ./ J0(1);
        dJdp_norm        = dJdp.'        .* p_val ./ J0(1);
        par_J_par_p_norm = par_J_par_p.' .* p_val ./ J0(1);
        par_x_star_par_p_norm = par_x_star_par_p.' .* p_val ./ x0(1:end-1,1);
        delta_p_norm = delta_p_change_activity  ./ p_val.';
    end
end

function [par_x_star_par_p_all_params, dJstar_dp_all_params, ...
            dJdp_all_params, par_J_par_p_all_params,...
            delta_p_change_activity_all_params] = local_sens_one_obj_all_param(x0, p, params, param_idxs, lambda, grad, hess, num_constr_nl, obj)
% construct matrix equation for a single J - slide 47 of lec 10 of MDO
% https://arc.aiaa.org/doi/10.2514/3.51191 - Sobieski 1982

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

    % left hand side matrix - independent of parameter
    matrix = local_sens_LHS_matrix(x0,p,hess,active,active_lin,active_lb,active_ub);
    
    % allocate outputs
    par_J_par_p_all_params = zeros(length(params),1);
    dJdp_all_params = zeros(length(params),1);
    dJstar_dp_all_params = zeros(length(params),1);
    par_x_star_par_p_all_params = zeros(length(params),length(x0)-1);
    num_constr_tot = length(lambda_nl) + length(lambda_lin) + length(lambda_lb) + length(lambda_ub);
    delta_p_change_activity_all_params = zeros(length(params),num_constr_tot);

    % sweep parameters
    parfor i = 1:length(params)
        param_name = params{i};
        param_idx = param_idxs(i);

        [par_J_par_p, dJdp, dJstar_dp, ...
         par_x_star_par_p, delta_p_change_activity] = local_sens_one_obj_one_param(matrix, ...
                                                        grad, obj, p, param_name, param_idx, x0, num_constr_nl, ...
                                                        lambda_nl, lambda_lin, lambda_lb, lambda_ub, ...
                                                        active, active_lin, active_lb, active_ub);

        % assign outputs
        par_J_par_p_all_params(i) = par_J_par_p;
        dJdp_all_params(i) = dJdp;
        dJstar_dp_all_params(i) = dJstar_dp;
        par_x_star_par_p_all_params(i,:) = par_x_star_par_p;
        delta_p_change_activity_all_params(i,:) = delta_p_change_activity;
    end
end

function [par_J_par_p, dJ_star_dp_lin, dJ_star_dp_quad, ...
    par_x_star_par_p, delta_p_change_activity] = local_sens_one_obj_one_param(matrix, grad, obj, ...
                                                        p, param_name, param_idx, x0, num_constr, ...
                                                        lambda_nl, lambda_lin, lambda_lb, lambda_ub, ...
                                                        active, active_lin, active_lb, active_ub)
        % right hand side vector - depends on parameter
        [vector, par_J_par_p, dJ_star_dp_lin] = local_sens_RHS_vector_and_dJdp(obj, p, param_name, param_idx, x0, num_constr, ...
                                                        lambda_nl, lambda_lin, lambda_lb, lambda_ub, ...
                                                        active, active_lin, active_lb, active_ub);
    
        % matrix solution
        sol = matrix \ vector;
        par_x_star_par_p = sol(1:length(x0)-1);
        par_lam_par_p = sol(length(x0):end);
    
        % total derivative
        par_J_par_x = grad;
        dJ_star_dp_quad = par_J_par_p + par_J_par_x.' * par_x_star_par_p;

        debug = false;
        if debug
            % comparision against global
            vars_global = load('global_sens.mat','X_LCOE','param_name','p','ratios');
            assert(strcmp(vars_global.param_name, param_name))
            dx = squeeze(vars_global.X_LCOE(1,4,:)  - vars_global.X_LCOE(1,3,:));
            dp = (vars_global.ratios(4)-vars_global.ratios(3)) * vars_global.p.(param_name);
            dxdp = dx(1:14) / dp;

            color_each_element(par_x_star_par_p)
            title('local dx*/dp')

            color_each_element(dxdp)
            title('global dx*/dp')

            pct_err = (par_x_star_par_p - dxdp)./dxdp*100;
            pct_err(abs(dxdp)<1e-4 & abs(par_x_star_par_p)<1e-4) = 0;
            color_each_element(pct_err)
            title('percent error')

            % breakdown of local into two terms
            M = inv(matrix);
            M1 = M(1:14,1:14);
            M2 = M(1:14,15:end);
            neg_c = vector(1:14);
            neg_d = vector(15:end);
            dxdp_term1 = M1 * neg_c;
            dxdp_term2 = M2 * neg_d;

            color_each_element(dxdp_term1)
            title('local dx*/dp term 1: -M_1 c')

            color_each_element(dxdp_term2)
            title('local dx*/dp term 2: -M_2 d')
        end

        % active constraints becoming inactive
        idx_nl = 1:sum(active);
        idx_lin = sum(active)+1 : sum(active)+sum(active_lin);
        idx_lb = sum(active)+sum(active_lin)+1 : sum(active)+sum(active_lin)+sum(active_lb);
        idx_ub = sum(active)+sum(active_lin)+sum(active_lb)+1 : sum(active)+sum(active_lin)+sum(active_lb)+sum(active_ub);
        delta_p_nl_becomes_inactive  = -lambda_nl(active) ./ par_lam_par_p(idx_nl);
        delta_p_lin_becomes_inactive = -lambda_lin(active_lin) ./ par_lam_par_p(idx_lin);
        delta_p_lb_becomes_inactive  = -lambda_lb(active_lb) ./ par_lam_par_p(idx_lb);
        delta_p_ub_becomes_inactive  = -lambda_ub(active_ub) ./ par_lam_par_p(idx_ub);

        % inactive constraints becoming active
        % fixme - this is zeroed for now
        delta_p_nl_becomes_active  = zeros(sum(~active),1);     %-g(~active) ./ (dg_dx(~active) * par_x_star_par_p)
        delta_p_lin_becomes_active = zeros(sum(~active_lin),1); %-g_lin(~active_lin) ./ (dg_lin_dx(~active_lin) * par_x_star_par_p)
        delta_p_lb_becomes_active  = zeros(sum(~active_lb),1);  %-g_lb(~active_lb) ./ (dg_lb_dx(~active_lb) * par_x_star_par_p)
        delta_p_ub_becomes_active  = zeros(sum(~active_ub),1);  %-g_ub(~active_ub) ./ (dg_ub_dx(~active_ub) * par_x_star_par_p)

        % combine constraint activity change for both active and inactive
        delta_p_nl_change_activity(active)  = delta_p_nl_becomes_inactive;
        delta_p_nl_change_activity(~active) = delta_p_nl_becomes_active;
        delta_p_lin_change_activity(active_lin)  = delta_p_lin_becomes_inactive;
        delta_p_lin_change_activity(~active_lin) = delta_p_lin_becomes_active;
        delta_p_lb_change_activity(active_lb)  = delta_p_lb_becomes_inactive;
        delta_p_lb_change_activity(~active_lb) = delta_p_lb_becomes_active;
        delta_p_ub_change_activity(active_ub)  = delta_p_ub_becomes_inactive;
        delta_p_ub_change_activity(~active_ub) = delta_p_ub_becomes_active;

        % combine all constraints: nonlin, lin, lb, ub
        delta_p_change_activity = [delta_p_nl_change_activity, delta_p_lin_change_activity, ...
                                   delta_p_lb_change_activity, delta_p_ub_change_activity];
end

function [vector, par_J_par_p, dJ_star_dp_lin] = local_sens_RHS_vector_and_dJdp(obj, p, param_name, param_idx, x0, num_constr, ...
                                                        lambda_nl, lambda_lin, lambda_lb, lambda_ub, ...
                                                        active, active_lin, active_lb, active_ub)


    p0 = p.(param_name);
    [par_J_par_p, par_g_par_p, ...
        par_par_J_par_x_par_p, ...
        par_par_g_par_x_par_p,...
        par_g_lin_par_p,...
        par_par_g_lin_par_x_par_p] = get_partials(obj, p, param_name, param_idx, x0, p0,num_constr);

    % assume that x bounds are parameter-independent
    par_g_lb_par_p = zeros(size(active_lb)); 
    par_g_ub_par_p = zeros(size(active_ub));

    % use all constraints for J* linear sensitivity - MDO book section 5.3.4, equation 5.35, page 175
    dJ_star_dp_lin = par_J_par_p ...
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
    num_constr_active_nl = sum(active);
    num_constr_active_tot = sum([active; active_lin; active_lb; active_ub]);

    % obtain B = dg/dx with finite difference for nonlinear constraints
    active_constraint_handle = @(X) get_constraints([X;1],p,active);
    B_nl = finite_difference_vector_x(active_constraint_handle,x0(1:end-1),num_constr_active_nl);

    % obtain B from constraint A_ineq matrix for linear constraints
    A_ineq_lin_part = lin_ineq_constraints(p);
    A_ineq_lin = zeros(size(A_ineq_lin_part,1), length(x0)-1);
    A_ineq_lin(:, 1:size(A_ineq_lin_part,2)) = A_ineq_lin_part;
    B_lin = A_ineq_lin(active_lin,:).';

    % B is +- 1 at for the active bounds
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
            par_par_g_lin_par_x_par_p] = get_partials(obj, p, param_name, param_idx, x0, p0,num_constr)

        %  compute par_y_par_p using finite difference
        y_fcn_handle = @(param_value) get_sim_outputs_and_derivs(obj,p,param_name,param_value,x0);
        if isscalar(p0)
            par_y_par_p = finite_difference_scalar_x(y_fcn_handle,p0);
        else
            par_y_par_p = finite_difference_scalar_x(y_fcn_handle,p0,param_idx);
        end
    
        % decompose par_y_par_p into its sections
        par_J_par_p = par_y_par_p(1,1);
        par_g_par_p = par_y_par_p(1,2 : 1+num_constr).';
        par_par_J_par_x_par_p = par_y_par_p(2:end, 1);
        par_par_g_par_x_par_p = par_y_par_p(2:end, 2 : 1+num_constr);

        % linear sensitivities (analytical)
        [~, ~, dAdp, dbdp] = lin_ineq_constraints(p, param_name);
        dAdp = dAdp(:,1:end-1); % remove material column
        par_g_lin_par_p = dAdp * x0(1:end-1) - dbdp;
        par_par_g_lin_par_x_par_p = dAdp.';
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
        x_before = x_1_vec(1:i-1);   % elements of x vector before the element being finite differenced
        x_after  = x_1_vec(i+1:end); % elements of x vector after  the element being finite differenced
        scalar_handle = @(x_i) fcn_handle([x_before; x_i; x_after]);
        deriv(i,:) = finite_difference_scalar_x(scalar_handle,x_1_vec(i));
    end
end

function deriv = finite_difference_scalar_x(f_handle,x_1,idx)
% fcn handle: a function expecting input scalar x and returning output
% matrix y.

% x is allowed to be nonscalar only if optional arg idx is passed, which
% indicates which element of x to use to find delta_x. The calculated
% delta_x is then applied to the entire x vector at once, so x still acts 
% like a scalar in terms of the derivative. This is different than 
% finite_difference_vector_x, which perturbs each element of x separately.

    if nargin <= 2
        assert(isscalar(x_1))
        idx = 1;
    end
    y_1 = f_handle(x_1);
    delta_x = max(1e-10, x_1(idx)*1e-4);
    x_2 = x_1 + delta_x;
    y_2 = f_handle(x_2);

    delta_y = y_2 - y_1;
    deriv = delta_y / delta_x;
    % deriv is size(y)
end

