function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         F_drag_f,F_drag_s,...
    phase_F_drag_f,phase_F_drag_s] = get_response_drag(H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                        gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,P_max,...
                                        drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                        X_max,control_type,control_solve_type,multibody,merge_bodies,...
                                        X_tol,phase_X_tol,F_lim_tol,X_lim_tol,...
                                        max_drag_iters_fixed_point,max_drag_iters_solver,...
                                        drag_convergence_plot_on,drag_fcn,...
                                        D_f,D_f_in,D_d,T_f_slam,T_s_slam)

    CTRL_MULT_MAX = 1e4;
    CTRL_ERR_REGULARIZATION = 1;

    % initial guess: 2m float amplitude, 0.5m spar amplitude
    X_f_guess = 2 * ones(size(w));
    phase_X_f_guess = zeros(size(w));
    X_s_guess = .5 * ones(size(w));
    phase_X_s_guess = zeros(size(w));
    if F_max == 0
        ctrl_mult_guess = zeros(size(w));
    else
        ctrl_mult_guess = ones(size(w));
    end
    phase_ctrl_mult_guess = zeros(size(w));

    % package inputs
    dynam_inputs = {H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,P_max,...
                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                    X_max,control_type,control_solve_type,multibody,merge_bodies,...
                    X_tol,phase_X_tol,F_lim_tol,X_lim_tol,...
                    drag_convergence_plot_on,drag_fcn,...
                    D_f,D_f_in,D_d,T_f_slam,T_s_slam};

    % first do fixed point iteration
    if max_drag_iters_fixed_point > 0
        [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         ctrl_mult,phase_ctrl_mult,...
         F_drag_f,F_drag_s,...
         phase_F_drag_f,phase_F_drag_s] = fixed_point_iteration(X_f_guess,X_s_guess,...
                                                    phase_X_f_guess,phase_X_s_guess,...
                                                    ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                    dynam_inputs{:}, ...
                                                    max_drag_iters_fixed_point,...
                                                    CTRL_ERR_REGULARIZATION);

        if strcmpi(control_solve_type,'solver')
            idx_use = ((ctrl_mult >= 0 & ctrl_mult < CTRL_MULT_MAX) & isfinite(mag_X_f)) | isnan(B_h_f);
        else
            idx_use = isfinite(mag_X_f) | isnan(B_h_f);
        end

        % update guesses
        [X_f_guess(idx_use),phase_X_f_guess(idx_use),...
         X_s_guess(idx_use),phase_X_s_guess(idx_use),...
         ctrl_mult_guess(idx_use),phase_ctrl_mult_guess(idx_use)] = deal(mag_X_f(idx_use),phase_X_f(idx_use),...
                                                       mag_X_s(idx_use),phase_X_s(idx_use),...
                                                       ctrl_mult(idx_use),phase_ctrl_mult(idx_use));

    end
    % then do nonlinear solver to finish
    if max_drag_iters_solver > 0
        
        [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         F_drag_f,F_drag_s,...
         phase_F_drag_f,phase_F_drag_s] = solver(X_f_guess,X_s_guess,...
                                        phase_X_f_guess,phase_X_s_guess,...
                                        ctrl_mult_guess,phase_ctrl_mult_guess,...
                                        dynam_inputs{:},max_drag_iters_solver);
    end

    % zero any negative power
    real_P(real_P<0) = 0;
end

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,...
     phase_F_drag_s] = solver(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                                                ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                                F_max,P_max,drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,control_solve_type,multibody,merge_bodies,...
                                                X_tol,phase_X_tol,F_lim_tol,X_lim_tol,...
                                                drag_convergence_plot_on,...
                                                drag_fcn,D_f,D_f_in,D_d,...
                                                T_f_slam,T_s_slam,max_drag_iters)

    % anonymous function to take collapsed row vector and turn into the nth sea state matrix
    sz = size(w);
    idx_not_nan = isfinite(B_h_f);
    [row,col] = find(idx_not_nan);
    if ~iscolumn(row) % ensure row and col are column vectors so accumarray works
        row = row.';
        col = col.';
    end
    N_ss_nz = nnz(idx_not_nan); % number of nonzero sea states (so not including NaNs)
    fill_nan_matrix = @(x) accumarray([row col], x, sz, [], NaN) ; 
    unflatten = @(x,n) fill_nan_matrix( x( (1:N_ss_nz) + (n-1)*N_ss_nz ) );

    if strcmpi(control_solve_type,'solver')
        num_solver_vars = 6;
    else
        num_solver_vars = 4;
    end

    % anonymous function to take sea state matrices and turn into collapsed row vector
    flatten = @(varargin) flatten_solver_vars(varargin, idx_not_nan, num_solver_vars);

    % prepare inputs for solver
    if num_solver_vars == 6
        x0 = flatten(X_f_guess, X_s_guess, phase_X_f_guess, phase_X_s_guess,...
                     ctrl_mult_guess, phase_ctrl_mult_guess);
    else
        x0 = flatten(X_f_guess, X_s_guess, phase_X_f_guess, phase_X_s_guess);
    end

    % fun_inner takes 1 input and returns the num_solver_vars+ outputs of dynamics_error_wrapper
    if num_solver_vars == 6
        fun_inner = @(x) dynamics_error_wrapper(unflatten(x,1),unflatten(x,2),...
                                    unflatten(x,3),unflatten(x,4),...
                                    unflatten(x,5),unflatten(x,6),...
                                    H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,P_max,...
                                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                    X_max,control_type,control_solve_type,multibody,merge_bodies,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam);
    else
        fun_inner = @(x) dynamics_error_wrapper(unflatten(x,1),unflatten(x,2),...
                                    unflatten(x,3),unflatten(x,4),...
                                    ctrl_mult_guess,phase_ctrl_mult_guess,...
                                    H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,P_max,...
                                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                    X_max,control_type,control_solve_type,multibody,merge_bodies,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam);
    end
    
    sparsity = repmat(eye(N_ss_nz),num_solver_vars);
    opts = optimoptions('fsolve','JacobPattern',sparsity,...
                                 'MaxIterations',max_drag_iters,...
                                 'Algorithm','trust-region',...
                                 ...%'StepTolerance',min([X_tol,phase_X_tol,F_lim_tol,X_lim_tol]),...
                                 'Display','off');
    if drag_convergence_plot_on
        opts.PlotFcn = {'optimplotx','optimplotfval'};
    end

    % solve
    [x_solved,~,~,~] = fsolve(@(x)fun_outer(x,fun_inner,flatten,num_solver_vars),x0,opts);

    % unpack
    X_f_solved = unflatten(x_solved,1);
    X_s_solved = unflatten(x_solved,2);
    phase_X_f_solved = unflatten(x_solved,3);
    phase_X_s_solved = unflatten(x_solved,4);
    if num_solver_vars == 6
        ctrl_mult_solved = unflatten(x_solved,5);
        phase_ctrl_mult_solved = unflatten(x_solved,6);
    else
        ctrl_mult_solved = ctrl_mult_guess;
        phase_ctrl_mult_solved = phase_ctrl_mult_guess;
    end

    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,...
     phase_F_drag_s] = dynamics_from_guess(X_f_solved, phase_X_f_solved, mag_v0_f, drag_const_f, ...
                                    X_s_solved, phase_X_s_solved, mag_v0_s, drag_const_s, ...
                                    B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,H,w,k_wvn,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                    control_type,multibody,merge_bodies,F_max,P_max,X_max,...
                                    control_solve_type,ctrl_mult_solved,phase_ctrl_mult_solved,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam, drag_convergence_plot_on);

end
function out_flat = fun_outer(x,fun_inner,flatten,num_solver_vars)
    Y = cell(1,num_solver_vars);
    [Y{:}] = fun_inner(x);
    out_flat = flatten(Y{:});
end

function out_flat = flatten_solver_vars(var_cell, idx_not_nan, num_solver_vars)
    chunks = cell(1,num_solver_vars);
    for i = 1:num_solver_vars
        chunks{i} = reshape(var_cell{i}(idx_not_nan),[],1);
    end
    out_flat = vertcat(chunks{:});
end

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     ctrl_mult_guess,...
     phase_ctrl_mult_guess,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,...
     phase_F_drag_s] = fixed_point_iteration(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                                                ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                                F_max,P_max,drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,control_solve_type,multibody,merge_bodies,...
                                                X_tol,phase_X_tol,F_lim_tol,X_lim_tol,...
                                                drag_convergence_plot_on,drag_fcn,...
                                                D_f, D_f_in, D_d, T_f_slam, ...
                                                T_s_slam, max_drag_iters, ctrl_err_regularization)

    converged = false;
    iters = 0;

    % loop until converged
    while ~converged
        if drag_convergence_plot_on
            % record guesses
            % 5,11 was found manually to be a problem sea state
            sea_state_plot = {5,11};
            X_f_guesses(iters+1) = X_f_guess(sea_state_plot{:});
            X_s_guesses(iters+1) = X_s_guess(sea_state_plot{:});
            phase_X_f_guesses(iters+1) = phase_X_f_guess(sea_state_plot{:});
            phase_X_s_guesses(iters+1) = phase_X_s_guess(sea_state_plot{:});
        end

        [X_f_err,X_s_err,...
         phase_X_f_err,phase_X_s_err,...
         force_lim_err,amp_lim_err,...
         mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         F_drag_f,F_drag_s,...
         phase_F_drag_f,phase_F_drag_s] = dynamics_error_wrapper(X_f_guess,X_s_guess,...
                                                       phase_X_f_guess,phase_X_s_guess,...
                                                       ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,P_max,...
                                                drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,control_solve_type,multibody,merge_bodies,...
                                                drag_fcn, D_f, D_f_in, D_d, ...
                                                T_f_slam, T_s_slam);

        X_err       = max( abs([X_f_err,X_s_err]),             [], 'all');
        phase_X_err = max( abs([phase_X_f_err,phase_X_s_err]), [], 'all');
        F_lim_err   = max( abs(force_lim_err),                 [], 'all');
        X_lim_err   = max( abs(amp_lim_err),                   [], 'all');

        if any(~isfinite(mag_X_f(:)) & isfinite(w(:)))
            warning('inf/nan in results but not in guess')
        end
        % new guesses
        X_f_guess = mag_X_f;
        X_s_guess = mag_X_s;
        phase_X_f_guess = phase_X_f;
        phase_X_s_guess = phase_X_s;
        if strcmpi(control_solve_type,'solver')
            % Heuristic update from prior solver-based flow: scale multiplier by
            % normalized constraint residual to improve next iterate.
            ctrl_mult_guess = ctrl_mult_guess ./ (force_lim_err + ctrl_err_regularization);
            phase_ctrl_mult_guess = zeros(size(ctrl_mult_guess));
        end

        % check convergence
        X_converged = X_err < X_tol;
        phase_X_converged = phase_X_err < phase_X_tol;
        if strcmpi(control_solve_type,'solver')
            F_lim_converged = F_lim_err < F_lim_tol;
            X_lim_converged = X_lim_err < X_lim_tol;
            converged = X_converged && phase_X_converged && ...
                        F_lim_converged && X_lim_converged;
        else
            converged = X_converged && phase_X_converged;
        end

        if all(isnan(X_f_guess))
            error('all nan')
        end

        % increment iterations and check max iters
        iters = iters + 1;
        if iters > max_drag_iters
            break
        end
    end

    if drag_convergence_plot_on % fixme: have guesses as output of dynamics and save it to vals and plot outside of the sim
        plot_drag_convergence(X_f_guesses, X_s_guesses, phase_X_f_guesses, phase_X_s_guesses, ...
                              iters, multibody, X_tol, phase_X_tol);
    end

end

function [X_f_err,X_s_err,...
          phase_X_f_err,phase_X_s_err,...
          force_lim_err,amp_lim_err,...
          mag_U,phase_U,...
          real_P,reactive_P,...
          mag_X_u,phase_X_u,...
          mag_X_f,phase_X_f,...
          mag_X_s,phase_X_s,...
          B_p,K_p,P_sat_ratio,...
          F_drag_f,F_drag_s,...
     phase_F_drag_f,phase_F_drag_s] = dynamics_error_wrapper(X_f_guess,X_s_guess,...
                                                        phase_X_f_guess,phase_X_s_guess,...
                                                        ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                                F_max,P_max,drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,control_solve_type,multibody,merge_bodies,drag_fcn,...
                                                D_f,D_f_in,D_d, T_f_slam, T_s_slam)
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,phase_F_drag_s,...
     force_lim_err,amp_lim_err] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                    X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                    B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,H,w,k_wvn,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                    control_type,multibody,merge_bodies,F_max,P_max,X_max,...
                                    control_solve_type,ctrl_mult_guess,phase_ctrl_mult_guess,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam, false); % closed loop dynamics; no polar plot during fsolve iterations
    % error
    X_f_err = X_f_guess - mag_X_f;
    X_s_err = X_s_guess - mag_X_s;
    phase_X_f_err = wrapToPi(phase_X_f_guess - phase_X_f); % wrapping error -pi to pi maintains continuity at zero,
    phase_X_s_err = wrapToPi(phase_X_s_guess - phase_X_s); % allowing for the solution to oscillate between +pi and -pi

    % Fixme: when response = Inf (closed loop unstable), this should be considered 
    % controller error when it's the controller's fault (system stabilizable)
    % and only zeroed when it's inevitable (system not stabilizable). 
    % For now, we zero these Infs in both cases.
    % Also handle finite precision errors that return nan.
    idx_unstable_or_fp_error = mag_X_f==Inf | (isnan(mag_X_f) & ~isnan(w));
    if any(idx_unstable_or_fp_error(:))
        warning('unstable')
    end
    idx_nan_power = isnan(real_P) & ~isnan(w);
    if any(idx_nan_power(:))
        warning('nan power for non-nan sea state')
    end
    X_f_err(idx_unstable_or_fp_error) = 0;
    X_s_err(idx_unstable_or_fp_error) = 0;
    phase_X_f_err(idx_unstable_or_fp_error) = 0;
    phase_X_s_err(idx_unstable_or_fp_error) = 0;
    force_lim_err(idx_unstable_or_fp_error) = 0;
    amp_lim_err(idx_unstable_or_fp_error) = 0;
end

function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p_sat,K_p_sat,...
         P_sat_ratio,...
     mag_F_drag_f,mag_F_drag_s,...
     phase_F_drag_f,phase_F_drag_s,...
         force_lim_err,...
         amp_lim_err] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                        X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                        B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,H,w,k_wvn,...
                                        gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                        control_type,multibody,merge_bodies,F_max,P_max,X_max,...
                                        control_solve_type,ctrl_mult_guess,phase_ctrl_mult_guess,...
                                        drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam, brute_force_plot_on)

    % linear system
    [control_evaluation_fcn,...
     real_G_u, imag_G_u, ...
     B_drag_f, B_drag_s, ...
     gamma_drag_f, gamma_drag_s] = linearize_around_guessed_response(X_f_guess, X_s_guess, ...
                                                                phase_X_f_guess, phase_X_s_guess, ...
                                                                mag_v0_f, mag_v0_s, H, w, k_wvn, ...
                                                                drag_const_f, drag_const_s,...
                                                                B_h_f, B_h_s, B_c, K_h_f, K_h_s, ...
                                                                m_c, m_f, m_s, gamma_f_mag, ...
                                                                gamma_f_phase, gamma_s_mag, ...
                                                                gamma_s_phase, multibody, merge_bodies,...
                                                                drag_fcn, D_f, D_f_in, D_d);

    
    % unsaturated optimal control gains
    [B_p,K_p] = controller(real_G_u, imag_G_u, w, control_type);
    stabilize_B = strcmpi(control_type,'reactive') || strcmpi(control_type,'damping');
    stabilize_K = strcmpi(control_type,'reactive');

    % unsaturated response (stabilized)
    [mag_U_unsat,P_unsat,...
     mag_X_u_unsat,...
     B_p_stabilized,K_p_stabilized,...
     mag_X_f_unsat,...
     mag_X_s_unsat] = control_evaluation_fcn(K_p,B_p,stabilize_B,stabilize_K);
    
    Z_th = 1 ./ complex(real_G_u, imag_G_u);

    if strcmpi(control_type,'none')
        % if no control, no control authority to enforce constraints
        mag_U = mag_U_unsat;
        real_P = P_unsat;
        mag_X_u = mag_X_u_unsat;
        mag_X_f = mag_X_f_unsat;
        mag_X_s = mag_X_s_unsat;
        B_p_sat = B_p_stabilized;
        K_p_sat = K_p_stabilized;

        % phases take a lot longer to compute, so only evaluate them for
        % the no control case - unsat phases aren't needed for the sat case
        [~,~,~,~,~,~,~,...
         phase_X_f,phase_X_s,...
         phase_U,reactive_P,phase_X_u] = control_evaluation_fcn(K_p,B_p,stabilize_B,stabilize_K);
        force_lim_err = zeros(size(w));
        amp_lim_err = zeros(size(w));
    elseif strcmpi(control_solve_type,'solver')
        [mag_U,real_P,...
         mag_X_u,mag_X_f,mag_X_s,...
         B_p_sat,K_p_sat,...
         force_lim_err, ...
         ~,...
         ~,~,...
         phase_X_f,phase_X_s,...
         amp_lim_err,...
         phase_U,reactive_P,phase_X_u] = response_and_ctrl_err_from_ctrl_guess(...
                                                    ctrl_mult_guess,...
                                                    phase_ctrl_mult_guess,...
                                                    F_max, P_max, mag_U_unsat,...
                                                    X_max, mag_X_f_unsat, ...
                                                    mag_X_s_unsat, ...
                                                    mag_X_u_unsat, ...
                                                    P_unsat, complex(B_p_stabilized, -K_p_stabilized ./ w), Z_th, w,...
                                                    control_evaluation_fcn,...
                                                    H,k_wvn,D_f,D_d,...
                                                    T_f_slam,T_s_slam,...
                                                    control_type);
    else
        % solve constrained optimal control via brute force
        [mag_U,phase_U,...
        real_P,reactive_P,...
        mag_X_u,phase_X_u,...
        mag_X_f,phase_X_f,...
        mag_X_s,phase_X_s,...
        B_p_sat,K_p_sat] = solve_brute_force_opt_control(Z_th, w, ...
                            mag_X_u_unsat, mag_X_f_unsat, mag_X_s_unsat, ...
                            X_f_guess, phase_X_f_guess, X_s_guess, phase_X_s_guess, ...
                            mag_U_unsat, P_unsat, B_p_stabilized, K_p_stabilized, ...
                             F_max, X_max, P_max, control_type, control_evaluation_fcn, ...
                             H, k_wvn, D_f, D_d, T_f_slam, T_s_slam, brute_force_plot_on && ~merge_bodies);
        force_lim_err = zeros(size(w));
        amp_lim_err = zeros(size(w));
    end

    P_sat_ratio = real_P ./ P_unsat;
    F_drag_f = -B_drag_f .* mag_X_f .* w .* exp(1i * (phase_X_f + pi/2)) + gamma_drag_f .* H/2;
    F_drag_s = -B_drag_s .* mag_X_s .* w .* exp(1i * (phase_X_s + pi/2)) + gamma_drag_s .* H/2;
    mag_F_drag_f = abs(F_drag_f);
    mag_F_drag_s = abs(F_drag_s);
    phase_F_drag_f = my_angle(F_drag_f);
    phase_F_drag_s = my_angle(F_drag_s);
end

function [opt_mag_U,opt_phase_U,...
    opt_real_P,opt_reactive_P,...
    opt_mag_X_u,opt_phase_X_u,...
    opt_mag_X_f,opt_phase_X_f,...
    opt_mag_X_s,opt_phase_X_s,...
    opt_B_p_sat,opt_K_p_sat] = solve_brute_force_opt_control(Z_th, w, ...
                            mag_X_u_unsat, mag_X_f_unsat, mag_X_s_unsat, ...
                            X_f_guess, phase_X_f_guess, X_s_guess, phase_X_s_guess, ...
                            mag_U_unsat, P_unsat, B_p_unsat, K_p_unsat, ...
                            F_max, X_max, P_max, control_type, ...
                            control_evaluation_fcn, H, k_wvn, D_f, D_d,...
                            T_f_slam,T_s_slam, brute_force_plot_on)
    
    Z_p_unsat = complex(B_p_unsat, - K_p_unsat ./ w);

    size_opt_ctrl_mesh = 11;  
    if rem(size_opt_ctrl_mesh,2)==0 % input is even
        size_mult_mesh  = size_opt_ctrl_mesh + 1; % this MUST be odd so ctrl_mult_guess includes 1!
        size_phase_mesh = size_opt_ctrl_mesh;    % this should ideally be even for symmetrical plots
    else
        size_mult_mesh  = size_opt_ctrl_mesh;
        size_phase_mesh = size_opt_ctrl_mesh + 1;
    end
    ctrl_mult_guess = logspace(-1,1,size_mult_mesh);
    if strcmp(control_type,'reactive')
        phase_ctrl_mult_guess = linspace(0,2*pi,size_phase_mesh+1);
        phase_ctrl_mult_guess = phase_ctrl_mult_guess(1:end-1); % 2pi=0 so don't repeat it
    else
        phase_ctrl_mult_guess = 0;
    end
    [MAG_GUESS,PHASE_GUESS] = meshgrid(ctrl_mult_guess,phase_ctrl_mult_guess);

    mag_U = nan([size(w), numel(MAG_GUESS)]);
    real_P = nan([size(w), numel(MAG_GUESS)]);
    mag_X_u = nan([size(w), numel(MAG_GUESS)]);
    mag_X_f = nan([size(w), numel(MAG_GUESS)]);
    mag_X_s = nan([size(w), numel(MAG_GUESS)]);
    B_p_sat = nan([size(w), numel(MAG_GUESS)]);
    K_p_sat = nan([size(w), numel(MAG_GUESS)]);
    constraint_err = nan([size(w), numel(MAG_GUESS)]);
    mag_ctrl_mult_stabilized = nan([size(w), numel(MAG_GUESS)]);
    phase_ctrl_mult_stabilized = nan([size(w), numel(MAG_GUESS)]);

    for i = 1:numel(MAG_GUESS)
        [mag_U(:,:,i),real_P(:,:,i),...
        mag_X_u(:,:,i),mag_X_f(:,:,i),mag_X_s(:,:,i),...
        B_p_sat(:,:,i),K_p_sat(:,:,i),...
        constraint_err(:,:,i),...
        mag_ctrl_mult_stabilized(:,:,i),...
        phase_ctrl_mult_stabilized(:,:,i)] = response_and_ctrl_err_from_ctrl_guess(MAG_GUESS(i)*ones(size(w)),...
                                                    PHASE_GUESS(i)*ones(size(w)),...
                                                    F_max, P_max, mag_U_unsat,...
                                                    X_max, mag_X_f_unsat, ...
                                                    mag_X_s_unsat, ...
                                                    mag_X_u_unsat, ...
                                                    P_unsat, Z_p_unsat, Z_th, w,...
                                                    control_evaluation_fcn,...
                                                    H,k_wvn,D_f,D_d,...
                                                    T_f_slam,T_s_slam,...
                                                    control_type,...
                                                    X_f_guess, phase_X_f_guess,...
                                                    X_s_guess, phase_X_s_guess);

    end

    % save unmodified power for polar plot before infeasible controllers are set to -1
    real_P_unmodified = real_P;

    % choose best controller for each sea state
    each_controller_feasible = constraint_err <= 0;
    ctrl_dim = 3;
    sea_state_dim =[1 2];
    each_sea_state_feasible = any(each_controller_feasible,ctrl_dim);
    real_P(each_sea_state_feasible & ~each_controller_feasible) = -1; % negative power for infeasible controllers in sea states where a feasbile controller exists, so they aren't chosen
    [opt_real_P,idx_opt] = max(real_P,[],ctrl_dim);

    sea_state_infeasible = ~each_sea_state_feasible & ~isnan(w);
    if any(sea_state_infeasible(:))
        
        constraint_err_pos_pwr = constraint_err;
        constraint_err_pos_pwr(real_P_unmodified<0) = Inf;
        
        [min_err,idx_least_err] = min(constraint_err_pos_pwr,[],ctrl_dim);
        idx_opt(sea_state_infeasible) = idx_least_err(sea_state_infeasible);
        % Use sub2ind for correct 3D indexing: idx_opt gives ctrl-dim index for
        % each (Hs,T) sea state, so plain linear indexing real_P_unmodified(idx_opt)
        % would silently access the wrong elements of the [Nh x NT x N_ctrl] array.
        sz = size(real_P_unmodified);
        [r,c] = ndgrid(1:sz(1), 1:sz(2));
        lin_idx = sub2ind(sz, r, c, idx_opt);
        opt_real_P = real_P_unmodified(lin_idx);

        debug_print = brute_force_plot_on;
        if debug_print
            warning(['no feasible controller found for %d/%d sea states, setting signals to '...
                'least error solution. You may want to adjust the ctrl_mult_guess to be wider'], ...
                sum(sea_state_infeasible(:)), sum(~isnan(w(:))))
            % diagnostic: identify which constraint(s) are not satisfied
            infeas_mask = ~each_sea_state_feasible;
            fprintf('  Brute force diagnostics for infeasible sea states:\n');
            fprintf('    min constraint_err across controllers: min=%.4g, max=%.4g\n', ...
                    min(min_err(infeas_mask)), max(min_err(infeas_mask)));
    
            % diagnostic: identify the number of sea states for which each
            % controller was able to satisfy all constraints
            num_ss_ctrl_feasible = reshape(squeeze(sum(each_controller_feasible,sea_state_dim)),size(MAG_GUESS));
            mag_phs_num_feasible = [MAG_GUESS(:) PHASE_GUESS(:) num_ss_ctrl_feasible(:)];
            fprintf('   Controllers and the number of sea states they are feasible: \n[mag,phase,number]=\n[%s]\n', ...
                formattedDisplayText(mag_phs_num_feasible))
    
            % evaluate per-constraint errors at the best guess for each infeasible sea state
            best_B_p = nan(size(w));
            best_mag_U_val = nan(size(w));
            best_mag_X_f_val = nan(size(w));
            best_mag_X_s_val = nan(size(w));
            best_mag_X_u_val = nan(size(w));
            best_real_P_val = nan(size(w));
            for idx_flat = find(infeas_mask(:))'
                linear_idx = idx_flat + (idx_least_err(idx_flat)-1)*numel(w);
                best_B_p(idx_flat) = B_p_sat(linear_idx);
                best_mag_U_val(idx_flat) = mag_U(linear_idx);
                best_mag_X_f_val(idx_flat) = mag_X_f(linear_idx);
                best_mag_X_s_val(idx_flat) = mag_X_s(linear_idx);
                best_mag_X_u_val(idx_flat) = mag_X_u(linear_idx);
                best_real_P_val(idx_flat) = real_P(linear_idx);
            end
            n_infeas = sum(infeas_mask(:));
            fprintf('    B_p<0 (neg damping):   %d/%d infeasible sea states\n', sum(best_B_p(infeas_mask) < 0), n_infeas);
            fprintf('    force exceeded F_max:  %d/%d (F_max=%.4g)\n', sum(best_mag_U_val(infeas_mask) > 4/pi*F_max), n_infeas, F_max);
            fprintf('    power exceeded P_max:  %d/%d (P_max=%.4g)\n', sum(best_real_P_val(infeas_mask) > P_max), n_infeas, P_max);
            fprintf('    float amp exceeded:    %d/%d (X_max(1)=%.4g)\n', sum(best_mag_X_f_val(infeas_mask) > X_max(1)), n_infeas, X_max(1));
            fprintf('    spar amp exceeded:     %d/%d (X_max(2)=%.4g)\n', sum(best_mag_X_s_val(infeas_mask) > X_max(2)), n_infeas, X_max(2));
            fprintf('    rel amp exceeded:      %d/%d (X_max(3)=%.4g)\n', sum(best_mag_X_u_val(infeas_mask) > X_max(3)), n_infeas, X_max(3));
        end
    end

    % polar plot of controller search grid (only in operational, non-mergedbodies case)
    if brute_force_plot_on
        make_ctrl_polar_plot(MAG_GUESS, PHASE_GUESS, real_P_unmodified, ...
                                constraint_err, idx_opt, ...
                                mag_ctrl_mult_stabilized, phase_ctrl_mult_stabilized);
    end

    % Re-evaluate the optimal controller once to obtain all phase outputs that were
    % skipped during the brute-force loop (phase_U, reactive_P, phase_X_u, phase_X_f,
    % phase_X_s are not needed to find the optimum and are only computed here).
    [opt_mag_U, ~, opt_mag_X_u, opt_mag_X_f, opt_mag_X_s, ...
     opt_B_p_sat, opt_K_p_sat, ~, ~, ~, ...
     opt_phase_X_f, opt_phase_X_s, ...
     ~, opt_phase_U, opt_reactive_P, opt_phase_X_u] = ...
        response_and_ctrl_err_from_ctrl_guess(MAG_GUESS(idx_opt), PHASE_GUESS(idx_opt), ...
                                              F_max, P_max, mag_U_unsat, X_max, ...
                                              mag_X_f_unsat, mag_X_s_unsat, mag_X_u_unsat, ...
                                              P_unsat, Z_p_unsat, Z_th, w, ...
                                              control_evaluation_fcn, H, k_wvn, D_f, D_d, ...
                                              T_f_slam, T_s_slam, control_type);

    % Warn if the final optimal response contains NaN for valid sea states.
    % This path is not covered by the warning in dynamics_error_wrapper because
    % that function is only called during the iterative drag-convergence loop.
    if any(isnan(opt_mag_X_f(~isnan(w))),'all')
        warning('NaN in opt_mag_X_f after brute-force re-evaluation for non-NaN sea states')
    end
    if any(isnan(opt_real_P(~isnan(w))),'all')
        warning('NaN in opt_real_P after brute-force re-evaluation for non-NaN sea states')
    end
end

function [mag_U,real_P,...
        mag_X_u,mag_X_f,mag_X_s,...
        B_stabilized,K_stabilized,...
        constr_viol_err,...
        mag_ctrl_mult_stabilized,...
        phase_ctrl_mult_stabilized,...
        phase_X_f,phase_X_s,...
        optimality_err,...
        phase_U,reactive_P,phase_X_u] = response_and_ctrl_err_from_ctrl_guess(...
                                                    ctrl_mult_guess,...
                                                    phase_ctrl_mult_guess,...
                                                    F_max, P_max, mag_U_unsat,...
                                                    X_max, mag_X_f_unsat, ...
                                                    mag_X_s_unsat, ...
                                                    mag_X_u_unsat, ...
                                                    real_P_unsat,Z_p, Z_th, w,...
                                                    control_evaluation_fcn,...
                                                    H,k,D_f,D_d,T_f,T_s,...
                                                    control_type,...
                                                    mag_X_f_guess, phase_X_f_guess,...
                                                    mag_X_s_guess, phase_X_s_guess)

    % saturated control gains from guess
    mult = ctrl_mult_guess .* exp(1i*phase_ctrl_mult_guess);
    Z_p_sat = Z_p .* mult;
    B_p_sat =       real(Z_p_sat);
    K_p_sat = -w .* imag(Z_p_sat);

    % saturated response (with B stabilization, and K stabilization for reactive control)
    stabilize_B = true;
    stabilize_K = strcmpi(control_type,'reactive');

    % Determine whether pre-computed X_f/X_s guesses were supplied (fast path).
    % When nargin > 21 the caller passes guessed mag/phase for X_f and X_s so
    % that the brute-force loop only needs to evaluate power, mag_U, and mag_X_u
    % from the controller (outputs 1-5 of evaluate_any_controller after the
    % B_l/K_l reorder).  The guesses avoid calling multibody_response for the
    % X_f/X_s rows and their phases in every loop iteration.
    use_X_guess = nargin > 21;

    if use_X_guess
        % Fast path: only evaluate power, mag_U, mag_X_u, B_stabilized, K_stabilized.
        % Outputs 1-5 of evaluate_any_controller after the B_l/K_l reorder are
        % [mag_U, real_P, mag_X_u, B_l, K_l], so nargout=5 skips X_f/X_s.
        [mag_U,real_P,...
         mag_X_u,...
         B_stabilized,K_stabilized] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);
        mag_X_f  = mag_X_f_guess;
        phase_X_f = phase_X_f_guess;
        mag_X_s  = mag_X_s_guess;
        phase_X_s = phase_X_s_guess;
    elseif nargout > 13
        % caller needs phase_U, reactive_P, phase_X_u — request all outputs
        [mag_U,real_P,...
         mag_X_u,...
         B_stabilized,K_stabilized,...
         mag_X_f,mag_X_s,...
         phase_X_f,phase_X_s,...
         phase_U,reactive_P,...
         phase_X_u] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);
    else
        % caller only needs magnitudes and phase_X_f/phase_X_s (for slamming constraint)
        [mag_U,real_P,...
         mag_X_u,...
         B_stabilized,K_stabilized,...
         mag_X_f,mag_X_s,...
         phase_X_f,phase_X_s] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);
    end
    Z_p_stabilized = complex(B_stabilized, -K_stabilized ./ w);

    stabilization_mult = Z_p_stabilized ./ Z_p_sat;
    ctrl_mult_stabilized = mult .* stabilization_mult;
    mag_ctrl_mult_stabilized = abs(ctrl_mult_stabilized);
    phase_ctrl_mult_stabilized = my_angle(ctrl_mult_stabilized);

    % control error: uses saturated and unsaturated response
    if nargout > 12
        [constr_viol_err,optimality_err] = control_errors_from_sat_results(ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                    F_max,P_max,mag_U_unsat,...
                                                    X_max(1), mag_X_f_unsat, mag_X_f, phase_X_f,...
                                                    X_max(2), mag_X_s_unsat, mag_X_s, phase_X_s,...
                                                    X_max(3), mag_X_u_unsat,...
                                                    real_P_unsat, Z_th, ...
                                                    H, k, D_f, D_d, T_f, T_s,...
                                                    mag_U, mag_X_u, B_p_sat, real_P);
    else
        constr_viol_err = control_errors_from_sat_results(ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                    F_max,P_max,mag_U_unsat,...
                                                    X_max(1), mag_X_f_unsat, mag_X_f, phase_X_f,...
                                                    X_max(2), mag_X_s_unsat, mag_X_s, phase_X_s,...
                                                    X_max(3), mag_X_u_unsat,...
                                                    real_P_unsat, Z_th, ...
                                                    H, k, D_f, D_d, T_f, T_s,...
                                                    mag_U, mag_X_u, B_p_sat, real_P);
    end
end

function [control_evaluation_fcn,...
          real_G_u, imag_G_u,...
          B_drag_f, B_drag_s,...
          gamma_drag_f, gamma_drag_s] = linearize_around_guessed_response(X_f_guess, X_s_guess, ...
                                                                    phase_X_f_guess, phase_X_s_guess, ...
                                                                    mag_v0_f, mag_v0_s, H, w, k_wvn, ...
                                                                    drag_const_f, drag_const_s,...
                                                                    B_h_f, B_h_s, B_c, K_f, K_s, ...
                                                                    m_c, m_f, m_s, gamma_f_lin_mag, ...
                                                                    gamma_f_lin_phase, gamma_s_lin_mag, ...
                                                                    gamma_s_lin_phase, multibody, merge_bodies,...
                                                                    drag_fcn, D_f, D_f_in, D_d)
    % drag coefficients linearized around guess response
    [B_drag_f, gamma_drag_f] = get_drag_dynamic_coeffs(X_f_guess, phase_X_f_guess, mag_v0_f, w, drag_const_f, drag_fcn, k_wvn, H, D_f/2, D_f_in/2);
    [B_drag_s, gamma_drag_s] = get_drag_dynamic_coeffs(X_s_guess, phase_X_s_guess, mag_v0_s, w, drag_const_s, drag_fcn, k_wvn, H, D_d/2, 0);

    % total linear system coefficients
    B_f = B_h_f + B_drag_f;
    B_s = B_h_s + B_drag_s;
    
    gamma_f_lin = gamma_f_lin_mag .* exp(1i * gamma_f_lin_phase);
    gamma_s_lin = gamma_s_lin_mag .* exp(1i * gamma_s_lin_phase);

    gamma_f = gamma_f_lin + gamma_drag_f;
    gamma_s = gamma_s_lin + gamma_drag_s;

    F_f = gamma_f .* H / 2;               % excitation force of wave
    F_s = gamma_s .* H / 2;

    F_f_mag   = abs(F_f);
    F_s_mag   = abs(F_s);
    F_f_phase = my_angle(F_f);
    F_s_phase = my_angle(F_s);

    % thevenin admittance of the controlled DOF as seen from load
    if multibody
        [real_G_u,imag_G_u,D_sys] = multibody_impedance(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w);
        
%         X_f = m_f.*w - K_f./w;
%         X_s = m_s.*w - K_s./w;
%         X_c = m_c.*w;
% 
%         Z_f = B_f + 1i*X_f;
%         Z_s = B_s + 1i*X_s;
%         Z_c = B_c + 1i*X_c;
% 
%         det_Z = Z_f .* Z_s - Z_c.^2;
%         sum_Z = Z_f + Z_s + 2 * Z_c;
%         Z_th = det_Z ./ sum_Z;
%         G_u = 1./Z_th;
% 
%         mag_G_u_squared = abs(G_u).^2;
%         real_G_u = real(G_u);
%         imag_G_u = imag(G_u);
    else
        % derived on page 58 of notebook
        resistance = B_f;
        reactance = m_f.*w - K_f./w;
        mag_G_u_squared = 1 ./ (resistance.^2 + reactance.^2);
        real_G_u = resistance .* mag_G_u_squared;
        imag_G_u = -reactance .* mag_G_u_squared;
        D_sys = [];
    end

    idx_not_stabilizable = real_G_u < 0; % controlled dof sees negative source 
                                         % impedance, so it's not only unstable
                                         % but also not stabilizable

    % Pre-compute all controller-independent helper terms for the fast path of
    % multibody_response.  D_sys (already computed above) is passed in to avoid
    % recomputation inside get_multibody_helper_terms.
    if multibody
        [h_t18,h_t9,h_t108,h_t121,h_t_110_118,h_t144,h_t109,h_t119,h_D_sys, ...
         h_t145,h_t103,h_t126,h_t130,h_t129,h_t131,h_t102,h_t124] = ...
            get_multibody_helper_terms(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w, ...
                                       F_f_mag,F_f_phase,F_s_mag,F_s_phase,D_sys);
        helper_terms = {h_t18,h_t9,h_t108,h_t121,h_t_110_118,h_t144,h_t109,h_t119,h_D_sys, ...
                        h_t145,h_t103,h_t126,h_t130,h_t129,h_t131,h_t102,h_t124};
    else
        helper_terms = {};
    end

    % function to evaluate any controller for fixed linear system
    control_evaluation_fcn = @(K_p,B_p,stabilize_B,stabilize_K) evaluate_any_controller(...
                                        B_c, B_f, B_s, K_f, K_s, ...
                                        m_c, m_f, m_s, w, K_p, B_p, ...
                                        F_f_mag, F_f_phase, F_s_mag, ...
                                        F_s_phase, multibody, merge_bodies,...
                                        idx_not_stabilizable, ...
                                        stabilize_B, stabilize_K, D_sys, helper_terms);
end

function [B_p,K_p] = controller(real_G_u, imag_G_u, w, control_type)

    mag_G_u_squared = real_G_u.^2 + imag_G_u.^2;

    % control - set powertrain coefficients
    if strcmpi(control_type,'reactive')
        K_p = -w .* imag_G_u ./ mag_G_u_squared;
        B_p = real_G_u ./ mag_G_u_squared;
    elseif strcmpi(control_type, 'damping')
        B_p = 1 ./ sqrt(mag_G_u_squared);
        K_p = 0;
    elseif strcmpi(control_type, 'none')
        B_p = 0;
        K_p = 0;
    end
    B_p = max(B_p,0); % prevent negative controller damping

end

function [B_drag_2, gamma_drag] = get_drag_dynamic_coeffs(X_guess, phase_X_guess, mag_v0, w, drag_const, drag_fcn, k_wvn, H, R, R_in)

    idx_inf = isinf(X_guess); % override unstable guesses to prevent extra nans
    X_guess(idx_inf) = 1;
    phase_X_guess(idx_inf) = 0;

    idx_neg = X_guess < 0;
    X_guess(idx_neg) = 0.01;

    mag_v = w .* X_guess;
    phase_v = phase_X_guess + pi/2;
    mag_v0_v_ratio = mag_v0 ./ mag_v;
    mag_v0_v_ratio(mag_v==0) = 0; % prevent divide by zero
    
    % integral method: force should equal above method when k is small
    % (long waves, large T). We don't expect old Bd=new Bd because they are
    % defined differently (how much of K vs gamma goes into Bd).
    r = mag_v0_v_ratio;
    kappa = k_wvn * R;
    theta = phase_v;
    alpha = R_in / R;
    [B_int_1,G_int_real_1,G_int_imag_1] = drag_fcn(r,theta,kappa);
    [B_int_2,G_int_real_2,G_int_imag_2] = drag_fcn(r,theta,alpha*kappa);
    % derived p127 of notebook 9 2/25/26
    B_integral_weighted      = B_int_1      - alpha.^2 .* B_int_2;
    G_integral_real_weighted = G_int_real_1 - alpha.^2 .* G_int_real_2;
    G_integral_imag_weighted = G_int_imag_1 - alpha.^2 .* G_int_imag_2;
    G_integral_weighted = complex(G_integral_real_weighted, G_integral_imag_weighted);

    % excitation only when mag_v = 0
    mag_v_or_v0 = mag_v;
    idx_v_0 = mag_v == 0;
    mag_v_or_v0(idx_v_0) = mag_v0(idx_v_0);

    % derived p141 of notebook 9 3/5/26
    mag_v_constant_term = drag_const / (pi*(1-alpha^2)) * 2 * mag_v_or_v0;
    B_drag_2   = mag_v_constant_term           .* B_integral_weighted;
    assert(all( B_drag_2(isfinite(B_drag_2)) >= 0 ),'negative damping from integral method')
    gamma_drag = mag_v_constant_term .* (mag_v0 ./ (H/2)) .* G_integral_weighted;

    plot_on = false;
    drag_debug = false; % set true to plot new integral drag vs old drag
    if drag_debug
        % the following code makes it plot only once solver has converged
        converged = any(strcmp({dbstack().name},'solver')   & [dbstack().line]==166);
        op_seas   = any(strcmp({dbstack().name},'dynamics') & [dbstack().line]==72);
        if converged && op_seas
            plot_on = true;
        end
    end
    if plot_on
        phase_v_prime = atan2( cos(phase_X_guess) - mag_v0_v_ratio, -sin(phase_X_guess)); % derived on p67 of notebook #6 (10/4/24)
    
        alpha_v = sqrt(1 + mag_v0_v_ratio.^2 - 2 * mag_v0_v_ratio .* cos(phase_X_guess)); % derived on p48 of notebook #5 (6/7/24)
        phi_alpha = phase_v_prime; % - phase_v;
        mag_cf = drag_const * alpha_v.^2 .* mag_v; % eq 52 in Water paper

        B_drag = -1*mag_cf .* cos(phi_alpha); % real part of c_f
    %     B_drag = max(B_drag,zeros(size(B_drag))); % prevent negative damping
        K_drag = - w .* mag_cf .* sin(phi_alpha); % -w times imag part of c_f

        plot_drag_integral_debug(B_drag, mag_v, phase_v, K_drag, X_guess,...
                                 phase_X_guess, B_drag_2, gamma_drag);
    end
end

function [mag_U,real_P,...
     mag_X_u,...
     B_l,K_l,...
     mag_X_f,mag_X_s,...
     phase_X_f,phase_X_s,...
     phase_U,reactive_P,...
     phase_X_u] = evaluate_any_controller(...
                            B_c, B_f, B_s, K_f, K_s, ...
                            m_c, m_f, m_s, w, K_l, B_l, ...
                            F_f_mag, F_f_phase, F_s_mag, F_s_phase, ...
                            multibody, merge_bodies, idx_not_stabilizable,...
                            stabilize_B,stabilize_K,D_sys,helper_terms)

    % D_sys is an optional precomputed complex system determinant produced by
    % multibody_impedance.  It is independent of the controller (B_p/K_p) and
    % is passed through to multibody_response to avoid redundant recomputation.
    if nargin < 21
        D_sys = [];
    end
    % helper_terms is an optional 1x17 cell of controller-independent
    % intermediate variables from get_multibody_helper_terms.  When non-empty
    % they are forwarded to multibody_response as args 17-33 (nargin==33 fast path).
    if nargin < 22
        helper_terms = {};
    end

    % check stability
    [idx_closed_loop_unstable,...
     recommended_increase_Bl,...
     recommended_increase_Kl] = check_cl_stability(B_c, B_f, B_s, K_f, K_s, ...
                                                    m_c, m_f, m_s, w, ...
                                                    K_l, B_l, ...
                                                    idx_not_stabilizable, ...
                                                    multibody, merge_bodies);
    % stabilize if specified
    B_l = B_l + stabilize_B * recommended_increase_Bl;
    K_l = K_l + stabilize_K * recommended_increase_Kl;
    
    % check that stabilizing worked
    if stabilize_K || stabilize_B
        [idx_closed_loop_unstable,...
             recommended_increase_Bl,...
             recommended_increase_Kl] = check_cl_stability(B_c, B_f, B_s, K_f, K_s, ...
                                                            m_c, m_f, m_s, w, ...
                                                            K_l, B_l, idx_not_stabilizable, ...
                                                            multibody, merge_bodies);

        need_more_stabilizing = stabilize_K * max(abs(recommended_increase_Kl),[],'all') ...
                              + stabilize_B * max(abs(recommended_increase_Bl),[],'all');


        if need_more_stabilizing~=0 % allow double stabilizing due to finite precision
            B_l = B_l + stabilize_B * recommended_increase_Bl;
            K_l = K_l + stabilize_K * recommended_increase_Kl;

            [idx_closed_loop_unstable,...
             recommended_increase_Bl,...
             recommended_increase_Kl] = check_cl_stability(B_c, B_f, B_s, K_f, K_s, ...
                                                            m_c, m_f, m_s, w, ...
                                                            K_l, B_l, idx_not_stabilizable, ...
                                                            multibody, merge_bodies);

            need_more_stabilizing = stabilize_K * max(abs(recommended_increase_Kl),[],'all') ...
                                  + stabilize_B * max(abs(recommended_increase_Bl),[],'all');
            if need_more_stabilizing~=0
                warning('Stabilizing did not work after 2 tries. This could be a finite precision issue.')
            end

        end

        idx_closed_loop_unstable(recommended_increase_Kl==0 & recommended_increase_Bl==0) = 0; % required for finite precision
    end

    % response
    if merge_bodies
        B_m = B_f + B_s + 2*B_c;
        m_m = m_f + m_s + 2*m_c;
        K_m = K_f + K_s;
        
        F_exc = F_f_mag .* exp(1i*F_f_phase) + F_s_mag .* exp(1i*F_s_phase);
        F_m_mag = abs(F_exc);
        F_m_phase = my_angle(F_exc);

        [mag_X_f,phase_X_f] = second_order_transfer_fcn(w, m_m, B_m, K_m, F_m_mag, F_m_phase);
        mag_X_s = mag_X_f;
        phase_X_s = phase_X_f;
        [mag_X_u,phase_X_u,real_P,reactive_P] = deal(0*phase_X_f);

        Z_f = B_f + 1i * (m_f.*w.^2 - K_f ./ w);
        Z_s = B_s + 1i * (m_s.*w.^2 - K_s ./ w);
        Z_c = B_c + 1i * (m_c.*w.^2);
        F_brake_f = 1i*w .* (Z_f + Z_c) .* mag_X_f .* exp(1i*phase_X_f) - F_f_mag .* exp(1i*F_f_phase);
        F_brake_s = 1i*w .* (Z_s + Z_c) .* mag_X_f .* exp(1i*phase_X_f) - F_s_mag .* exp(1i*F_s_phase);
        mag_U_f = abs(F_brake_f);
        phase_U_f = my_angle(F_brake_f);
        mag_U_s = abs(F_brake_s);
        phase_U_s = my_angle(F_brake_s);
        % fixme it should be mag_U = mag_U_s and phase_U = pi+phase_U_s,
        % take the max and use float phase for now
        mag_U = max(mag_U_f,mag_U_s);
        phase_U = phase_U_f;

    else
        if multibody
            if nargout > 9
                [mag_U,real_P,...
                 mag_X_u,mag_X_f,mag_X_s,...
                 phase_X_f,phase_X_s,...
                 phase_U,reactive_P,...
                 phase_X_u] = multibody_response(B_c, B_f, B_s, K_f, K_s, ...
                                                 m_c, m_f, m_s, w, ...
                                                 K_l, B_l, ...
                                                 F_f_mag, F_f_phase, ...
                                                 F_s_mag, F_s_phase, D_sys, helper_terms{:});
            elseif nargout > 7
                [mag_U,real_P,...
                 mag_X_u,mag_X_f,mag_X_s,...
                 phase_X_f,phase_X_s] = multibody_response(B_c, B_f, B_s, K_f, K_s, ...
                                                           m_c, m_f, m_s, w, ...
                                                           K_l, B_l, ...
                                                           F_f_mag, F_f_phase, ...
                                                           F_s_mag, F_s_phase, D_sys, helper_terms{:});
            elseif nargout > 5
                % caller needs mag_X_f/mag_X_s (positions 6/7) but not phases
                [mag_U,real_P,...
                 mag_X_u,mag_X_f,mag_X_s] = multibody_response(B_c, B_f, B_s, K_f, K_s, ...
                                                               m_c, m_f, m_s, w, ...
                                                               K_l, B_l, ...
                                                               F_f_mag, F_f_phase, ...
                                                               F_s_mag, F_s_phase, D_sys, helper_terms{:});
            else
                % fast path: caller only needs mag_U, real_P, mag_X_u, B_l, K_l (pos 1-5);
                % B_l/K_l come from the stabilization above, not from multibody_response.
                [mag_U,real_P,mag_X_u] = multibody_response(B_c, B_f, B_s, K_f, K_s, ...
                                                             m_c, m_f, m_s, w, ...
                                                             K_l, B_l, ...
                                                             F_f_mag, F_f_phase, ...
                                                             F_s_mag, F_s_phase, D_sys, helper_terms{:});
            end
        else
            b_sat = B_f + B_l;
            k_sat = K_f + K_l;
    
            [mag_X_f,phase_X_f] = second_order_transfer_fcn(w, m_f, b_sat, k_sat, F_f_mag, F_f_phase);
            mag_X_u = mag_X_f;   phase_X_u = phase_X_f;
            mag_X_s = 0*mag_X_f; phase_X_s = 0*phase_X_f;
            F_ptrain_over_x = sqrt( (B_l .* w).^2 + (K_l).^2 );
            mag_U = F_ptrain_over_x .* mag_X_f;
            
            phase_Z_u = atan2(-K_l./w, B_l);    % phase of control impedance
            phase_V_u = pi/2 + phase_X_u;               % phase of control velocity
            phase_U = phase_V_u + phase_Z_u;            % phase of control force
    
            real_P = 1/2 * B_l .* w.^2 .* mag_X_u.^2; % this is correct even if X and U are out of phase
            check_P = 1/2 * w .* mag_X_u .* mag_U .* cos(phase_U - phase_V_u); % so is this, they match
            reactive_P = zeros(size(real_P)); % fixme value is incorrect but output not used currently
        end
    end

    % set values to Inf when closed loop unstable
    mag_U(idx_closed_loop_unstable)      = Inf;
    real_P(idx_closed_loop_unstable)     = Inf;
    mag_X_u(idx_closed_loop_unstable)    = Inf;
    % mag_X_f/mag_X_s are at output positions 6/7; always computed for
    % non-multibody/merge_bodies paths; for multibody they require nargout>5.
    if nargout > 5 || ~multibody || merge_bodies
        mag_X_f(idx_closed_loop_unstable) = Inf;
        mag_X_s(idx_closed_loop_unstable) = Inf;
    end
    if nargout > 7 || ~multibody || merge_bodies
        phase_X_f(idx_closed_loop_unstable)  = Inf;
        phase_X_s(idx_closed_loop_unstable)  = Inf;
    end
    if nargout > 9 || ~multibody || merge_bodies
        phase_U(idx_closed_loop_unstable)    = Inf;
        reactive_P(idx_closed_loop_unstable) = Inf;
        phase_X_u(idx_closed_loop_unstable)  = Inf;
    end
end

function [idx_closed_loop_unstable,...
          recommended_increase_Bl,...
          recommended_increase_Kl] = check_cl_stability(B_c, B_f, B_s, K_f, K_s, ...
                            m_c, m_f, m_s, w, K_l_sat, B_l_sat, ...
                            idx_not_stabilizable, multibody, merge_bodies)
    % Check if the closed loop system is stable in all stabilizable sea states.
    % Must be checked individually for every sea state because eig doesn't 
    % accept arrays of more than 2 dimensions. This check must consider all
    % modes of the system (2 for multibody), not just the controlled mode.
    % Only checking closed loop means open loop unstable is allowed (intentionally).

    % The unstable sea states (idx_closed_loop_unstable) are then the union 
    % of unstable stabilizable sea states (idx_ctrl_makes_closed_loop_unstable)
    % and unstabilizable sea states (idx_not_stabilizable). The latter is 
    % from a previous check that Re(G_u)>0. idx_not_stabilizable = false 
    % means that even if the open loop system is unstable, there exists a 
    % controller that can stabilize it (ie, the unstable mode is controllable).

    recommended_increase_Bl = zeros(size(B_f));
    recommended_increase_Kl = zeros(size(B_f));

    Z_l = complex(B_l_sat, -K_l_sat ./ w);
    Z_p = Z_l; % fixme apply cascase here
    % uncomment and vectorize these once PTO dynamics are added
%     if multibody
%         cascade_VI_to_FXdot = eye(2);
%         T = [1, -1];
%     else
%         cascade_VI_to_FXdot = 1;
%         T = 1;
%     end
%     num = [1 0] * cascade_VI_to_FXdot * [Z_l(idx_H,idx_T); 1];
%     den = [0 1] * cascade_VI_to_FXdot * [Z_l(idx_H,idx_T); 1];
%     Z_pto = num ./ den;
%     Z_p = T.' * Z_pto * T; % transform PTO impedance onto body DOFs
    B_p =  real(Z_p);
    K_p = -imag(Z_p).*w;

    if merge_bodies
        ol_stable = (m_f + 2*m_c + m_s > 0) & (B_f + 2*B_c + B_s > 0) & (K_f + K_s > 0);
        cl_stable = ol_stable;
    else
        if multibody
            ol_stable = check_posdef_three_2x2s(m_f, m_c, m_c, m_s, ...
                                                B_f, B_c, B_c, B_s,...
                                                K_f, 0,   0,   K_s);
            cl_stable = check_posdef_three_2x2s(m_f,     m_c,     m_c,     m_s, ...
                                                B_f+B_p, B_c-B_p, B_c-B_p, B_s+B_p,...
                                                K_f+K_p,    -K_p,    -K_p, K_s+K_p);
        else 
            ol_stable = (m_f > 0) & (B_f     > 0) & (K_f     > 0);
            cl_stable = (m_f > 0) & (B_f+B_p > 0) & (K_f+K_p > 0);
        end
    end
    ol_unstable = ~ol_stable & ~isnan(w);
    if any(ol_unstable)
        % open loop matrices should all be stable, otherwise
        % indicates an error with the coefficients
        [idx_H,idx_T] = find(ol_unstable);
        if merge_bodies
            det_m_ol = m_f + 2*m_c + m_s;
            det_b_ol = B_f + 2*B_c + B_s;
            det_k_ol = K_f + K_s;
        else
            if multibody
                det_m_ol = m_f.*m_s - m_c.^2;
                det_b_ol = B_f.*B_s - B_c.^2;
                det_k_ol = K_f.*K_s;
            else
                det_m_ol = m_f;
                det_b_ol = B_f;
                det_k_ol = K_f;
            end
        end
        warning(['Open loop dynamics are unstable for idx_H=[%s], ' ...
               'idx_T=[%s]. det(M)=[%s], det(B)=[%s], det(K)=[%s].'], ...
            num2str(idx_H.'), ...
            num2str(idx_T.'), ...
            num2str(det_m_ol(ol_unstable).','%.2g '), ...
            num2str(det_b_ol(ol_unstable).','%.2g '), ...
            num2str(det_k_ol(ol_unstable).','%.2g '));
    end

    idx_ctrl_makes_closed_loop_unstable = ~cl_stable & ~isnan(w) & ~idx_not_stabilizable;

    if any(idx_ctrl_makes_closed_loop_unstable(:)) && ~merge_bodies
        % compute required change in control to create
        % stability - see notebook 9 p117 2/21/26
        if multibody
            det_B = (B_f+B_p).*(B_s+B_p) - (B_c-B_p).^2;
            sum_B = B_f + B_s + 2*B_c;
            denom_B = sum_B;

            det_K = (K_f+K_p).*(K_s+K_p) - (-K_p).^2;
            sum_K = K_f + K_s;
            denom_K = sum_K;
        else
            det_B = B_f + B_p;
            denom_B = 1;
            
            det_K = K_f + K_p;
            denom_K = 1;
        end

        rec_incr_Bp = abs(det_B) ./ denom_B * 1.01;
        %rec_incr_Bp(det_B==0) = eps;
        rec_incr_Bl = rec_incr_Bp; % fixme this should use cascade matrix

        rec_incr_Kp = abs(det_K) ./ denom_K * 1.01;
        %rec_incr_Kp(det_K==0) = eps;
        rec_incr_Kl = rec_incr_Kp; % fixme this should use cascade matrix

        idx_change_B = idx_ctrl_makes_closed_loop_unstable & det_B <= 0;
        idx_change_K = idx_ctrl_makes_closed_loop_unstable & det_K <= 0;
        recommended_increase_Bl(idx_change_B) = rec_incr_Bl(idx_change_B);
        recommended_increase_Kl(idx_change_K) = rec_incr_Kl(idx_change_K);
    end
                    
    idx_closed_loop_unstable = idx_not_stabilizable | idx_ctrl_makes_closed_loop_unstable;
end

function all_posdef = check_posdef_three_2x2s(a11,a12,a21,a22, b11,b12,b21,b22, c11,c12,c21,c22)
    all_posdef = ...
            a11 > 0 & (a11.*a22 - a21.*a12) >= 0 & ...
            b11 > 0 & (b11.*b22 - b21.*b12) >= 0 & ...
            c11 > 0 & (c11.*c22 - c21.*c12) >= 0;
end

function [constr_viol_err,optimality_err] = control_errors_from_sat_results(ctrl_mult_mag,...
                                                   ctrl_mult_phase, F_max, P_max, mag_U_unsat,...
                                                   X_f_upper_limit_static, mag_X_f_unsat, mag_X_f, phase_X_f,...
                                                   X_s_upper_limit_static, mag_X_s_unsat, mag_X_s, phase_X_s,...
                                                   X_u_upper_limit_static, mag_X_u_unsat,...
                                                   P_unsat, Z_th,...
                                                   H, k, D_f, D_d, T_f, T_s,...
                                                   mag_U, mag_X_u, B_p_sat, P_sat)

    use_amp_sat = ~(isinf(X_f_upper_limit_static) && isinf(X_s_upper_limit_static) && isinf(X_u_upper_limit_static));

    if use_amp_sat
        [X_f_lower_limit_dynamic,...
         X_f_upper_limit_dynamic,...
         idx_f_imag] = get_slamming_min_max(H/2, k, D_f, phase_X_f, T_f);
        X_f_lower_limit_dynamic(idx_f_imag) = -abs(imag(X_f_lower_limit_dynamic(idx_f_imag)));
        X_f_upper_limit_dynamic(idx_f_imag) = -abs(imag(X_f_upper_limit_dynamic(idx_f_imag)));

        [X_s_lower_limit_dynamic,...
         X_s_upper_limit_dynamic,...
         idx_s_imag] = get_slamming_min_max(H/2, k, D_d, phase_X_s, T_s);
        X_s_upper_limit_dynamic(idx_s_imag) = -abs(imag(X_s_upper_limit_dynamic(idx_s_imag)));
        X_s_lower_limit_dynamic(idx_s_imag) = -abs(imag(X_s_lower_limit_dynamic(idx_s_imag)));

        X_f_upper_limit = min(X_f_upper_limit_dynamic, X_f_upper_limit_static);
        X_s_upper_limit = min(X_s_upper_limit_dynamic, X_s_upper_limit_static);
        X_u_upper_limit = X_u_upper_limit_static;

        X_f_lower_limit = X_f_lower_limit_dynamic;
        X_s_lower_limit = X_s_lower_limit_dynamic;
    else
        X_f_upper_limit = Inf(size(mag_X_f));
        X_s_upper_limit = Inf(size(mag_X_s));
        X_u_upper_limit = Inf(size(mag_X_u));
        X_f_lower_limit = -Inf(size(mag_X_f));
        X_s_lower_limit = -Inf(size(mag_X_s));
    end

    % what the force, amp, and power would be if the force-sat, amp-sat,
    % and power-sat solutions applied - notebook p106 2/2/25
    force_saturated_U = min(4/pi * F_max,      mag_U_unsat);

    amp_up_saturated_X_f = min(X_f_upper_limit, mag_X_f_unsat);
    amp_dn_saturated_X_f = max(X_f_lower_limit, mag_X_f_unsat);

    amp_up_saturated_X_s = min(X_s_upper_limit, mag_X_s_unsat);
    amp_dn_saturated_X_s = max(X_s_lower_limit, mag_X_s_unsat);

    amp_saturated_X_u = min(X_u_upper_limit, mag_X_u_unsat);

    P_min = 0;
    power_up_saturated_P = min(P_max, P_unsat);
    power_dn_saturated_P = max(P_min, P_unsat);

    % indices where sat and unsat solutions violate which constraint
    idx_force_viol_unsat    = mag_U_unsat   > force_saturated_U;
    idx_amp_f_viol_unsat_up = mag_X_f_unsat > amp_up_saturated_X_f;
    idx_amp_f_viol_unsat_dn = mag_X_f_unsat < amp_dn_saturated_X_f;
    idx_amp_s_viol_unsat_up = mag_X_s_unsat > amp_up_saturated_X_s;
    idx_amp_s_viol_unsat_dn = mag_X_s_unsat < amp_dn_saturated_X_s;
    idx_amp_u_viol_unsat    = mag_X_u_unsat > amp_saturated_X_u;
    idx_power_viol_unsat_up = P_unsat       > P_max;
    idx_power_viol_unsat_dn = P_unsat       < P_min;
    
    num_constr_viol_unsat = idx_force_viol_unsat    + idx_amp_f_viol_unsat_up ...
                          + idx_amp_f_viol_unsat_dn + idx_amp_s_viol_unsat_up ...
                          + idx_amp_s_viol_unsat_dn + idx_amp_u_viol_unsat ...
                          + idx_power_viol_unsat_up + idx_power_viol_unsat_dn;
    all_ok_unsat       = num_constr_viol_unsat==0 & ~isnan(P_unsat);
    only_force_viol_unsat    = idx_force_viol_unsat    & num_constr_viol_unsat==1;
    only_amp_f_viol_unsat_up = idx_amp_f_viol_unsat_up & num_constr_viol_unsat==1;
    only_amp_f_viol_unsat_dn = idx_amp_f_viol_unsat_dn & num_constr_viol_unsat==1;
    only_amp_s_viol_unsat_up = idx_amp_s_viol_unsat_up & num_constr_viol_unsat==1;
    only_amp_s_viol_unsat_dn = idx_amp_s_viol_unsat_dn & num_constr_viol_unsat==1;
    only_amp_u_viol_unsat    = idx_amp_u_viol_unsat    & num_constr_viol_unsat==1;
    only_power_viol_unsat_up = idx_power_viol_unsat_up & num_constr_viol_unsat==1;
    only_power_viol_unsat_dn = idx_power_viol_unsat_dn & num_constr_viol_unsat==1;
    mult_const_viol_unsat = num_constr_viol_unsat >= 2;

    % indices where each solution applies. commented forumlas attempt to 
    % incorporate both but should not actually be based purely on whether
    % the saturated and/or unsaturated solution violates, it should also
    % depend on alpha and the limits (see notebook p142-144 9/15/25).
    idx_force_sat_applies    = all_ok_unsat | only_force_viol_unsat;
    idx_amp_f_sat_applies_up = all_ok_unsat | only_amp_f_viol_unsat_up;
    idx_amp_f_sat_applies_dn = all_ok_unsat | only_amp_f_viol_unsat_dn;
    idx_amp_s_sat_applies_up = all_ok_unsat | only_amp_s_viol_unsat_up;
    idx_amp_s_sat_applies_dn = all_ok_unsat | only_amp_s_viol_unsat_dn;
    idx_amp_u_sat_applies    = all_ok_unsat | only_amp_u_viol_unsat;
    idx_power_sat_applies_up = all_ok_unsat | only_power_viol_unsat_up;
    idx_power_sat_applies_dn = all_ok_unsat | only_power_viol_unsat_dn;

    % errors from each solution: >0 when constr is violated, and <0
    % when constr is not tight. 
    if isinf(F_max)
        F_err_from_force_sat = zeros(size(mag_X_f));
    else
        F_err_from_force_sat = mag_U   ./ force_saturated_U - 1;
        idx_U_0 = force_saturated_U == 0;
        F_err_from_force_sat(idx_U_0) = mag_U(idx_U_0) / 1e6;
        F_err_from_force_sat(F_err_from_force_sat < -1) = -.99;
    end
    X_err_from_amp_f_sat_up = mag_X_f ./ amp_up_saturated_X_f - 1;
    idx_amp_f_0 = amp_up_saturated_X_f == 0;
    X_err_from_amp_f_sat_up(idx_amp_f_0) = mag_X_f(idx_amp_f_0);
    X_err_from_amp_f_sat_up(X_err_from_amp_f_sat_up < -1) = -.99;

    X_err_from_amp_f_sat_dn = amp_dn_saturated_X_f ./ mag_X_f  - 1;
    idx_X_f_0 = mag_X_f == 0;
    X_err_from_amp_f_sat_dn(idx_X_f_0) = -amp_dn_saturated_X_f(idx_X_f_0);
    X_err_from_amp_f_sat_dn(X_err_from_amp_f_sat_dn < -1) = -.99;

    X_err_from_amp_s_sat_up = mag_X_s ./ amp_up_saturated_X_s - 1;
    idx_amp_s_0 = amp_up_saturated_X_s == 0;
    X_err_from_amp_s_sat_up(idx_amp_s_0) = mag_X_s(idx_amp_s_0);
    X_err_from_amp_s_sat_up(X_err_from_amp_s_sat_up < -1) = -.99;

    X_err_from_amp_s_sat_dn = amp_dn_saturated_X_s ./ mag_X_s  - 1;
    idx_X_s_0 = mag_X_s == 0;
    X_err_from_amp_s_sat_dn(idx_X_s_0) = -amp_dn_saturated_X_s(idx_X_s_0);
    X_err_from_amp_s_sat_dn(X_err_from_amp_s_sat_dn < -1) = -.99;

    X_err_from_amp_u_sat = mag_X_u ./ amp_saturated_X_u - 1;
    idx_amp_u_0 = amp_saturated_X_u == 0;
    X_err_from_amp_u_sat(idx_amp_u_0) = mag_X_u(idx_amp_u_0);
    X_err_from_amp_u_sat(X_err_from_amp_u_sat < -1) = -.99;

    if isinf(P_max)
        P_err_from_power_sat_up = zeros(size(mag_X_f));
    else
        P_err_from_power_sat_up = P_sat ./ power_up_saturated_P - 1;
        idx_P_up_0 = power_up_saturated_P == 0;
        P_err_from_power_sat_up(idx_P_up_0) = P_sat(idx_P_up_0) / 1e6;
        P_err_from_power_sat_up(P_err_from_power_sat_up < -1) = -.99;
    end

    P_err_from_power_sat_dn = power_dn_saturated_P ./ P_sat - 1;
    idx_P_dn_0 = P_sat == 0;
    P_err_from_power_sat_dn(idx_P_dn_0) = -power_dn_saturated_P(idx_P_dn_0) / 1e6;
    P_err_from_power_sat_dn(P_err_from_power_sat_dn < -1) = -.99;

    if ~use_amp_sat
        X_err_from_amp_f_sat_up = zeros(size(mag_X_f));
        X_err_from_amp_f_sat_dn = zeros(size(mag_X_f));
        X_err_from_amp_s_sat_up = zeros(size(mag_X_s));
        X_err_from_amp_s_sat_dn = zeros(size(mag_X_s));
        X_err_from_amp_u_sat    = zeros(size(mag_X_u));
    end

    constr_viol_err = Inf(size(mag_X_f));
    

    % at sea states where one of the saturated solutions necessarily applies, 
    % set constraint violation error as deviation from that solution,
    % and optimality error as deviation from the optimal ctrl_mult_phase.
    % 
    constr_viol_err(idx_force_sat_applies)    = F_err_from_force_sat(idx_force_sat_applies);
    constr_viol_err(idx_amp_f_sat_applies_up) = X_err_from_amp_f_sat_up(idx_amp_f_sat_applies_up);
    constr_viol_err(idx_amp_f_sat_applies_dn) = X_err_from_amp_f_sat_dn(idx_amp_f_sat_applies_dn);
    constr_viol_err(idx_amp_s_sat_applies_up) = X_err_from_amp_s_sat_up(idx_amp_s_sat_applies_up);
    constr_viol_err(idx_amp_s_sat_applies_dn) = X_err_from_amp_s_sat_dn(idx_amp_s_sat_applies_dn);
    constr_viol_err(idx_amp_u_sat_applies)    = X_err_from_amp_u_sat(idx_amp_u_sat_applies);
    constr_viol_err(idx_power_sat_applies_up) = P_err_from_power_sat_up(idx_power_sat_applies_up);
    constr_viol_err(idx_power_sat_applies_dn) = P_err_from_power_sat_dn(idx_power_sat_applies_dn);
    
    % at sea states where none of the solutions necessarily apply (the true 
    % solution could be any combo of one or two of the saturated solutions but unclear which 1-2),
    % set one error as the avg constraint violation (to penalize exceeding the constraint
    % without penalizing going under the constraint, so does not maximize power)
    % and the other error as the deviation from the max-power point when no
    % constraints are active or zero when a constraint is active (so that
    % when in the feasible region, it walks toward the max-power point
    % until it hits a constraint). To prevent it from starting infeasible
    % and terminating at the nearest constraint without ever going to the 
    % interior power-maximizing part, set optimality error=10 when any constraint is
    % violated. This will make it take larger steps when approaching from 
    % the outside so hopefully it overshoots and re-approaches from the inside.

    if isinf(P_max)
        constr_violation = (max(F_err_from_force_sat,   0) + max(X_err_from_amp_f_sat_up,0) + ...
                            max(X_err_from_amp_f_sat_dn,0) + max(X_err_from_amp_s_sat_up,0) + ...
                            max(X_err_from_amp_s_sat_dn,0) + max(X_err_from_amp_u_sat,   0) + ...
                            max(P_err_from_power_sat_dn))  / 7;
    else
        constr_violation = (max(F_err_from_force_sat,   0) + max(X_err_from_amp_f_sat_up,0) + ...
                            max(X_err_from_amp_f_sat_dn,0) + max(X_err_from_amp_s_sat_up,0) + ...
                            max(X_err_from_amp_s_sat_dn,0) + max(X_err_from_amp_u_sat,   0) + ...
                            max(P_err_from_power_sat_dn,0) + max(P_err_from_power_sat_up,0)   )  / 8;
    end
    constr_viol_err(mult_const_viol_unsat) = constr_violation(mult_const_viol_unsat);
    
    % prevent 0/0=NaN when limit of zero is set
    if F_max==0
        constr_viol_err(mag_U==0) = 0;
        constr_viol_err(mag_U~=0) = mag_U(mag_U~=0) ./ mag_U_unsat(mag_U~=0);
    else
        % add penalty for B_p<0 (negative power)
        tol = 0;
        B_p_violation = max(-B_p_sat+tol,0);
        constr_viol_err = constr_viol_err + B_p_violation;
    end

    if any(~isfinite(constr_viol_err(~isnan(H))),'all')
        warning('error is non finite for non-nan sea state')
    end
    if nargout > 1
        optimality_err  = Inf(size(mag_X_f));

        idx_soln_applies = idx_force_sat_applies    | idx_amp_f_sat_applies_up |...
                           idx_amp_f_sat_applies_dn | idx_amp_s_sat_applies_up |...
                           idx_amp_s_sat_applies_dn | idx_amp_u_sat_applies    | ...
                           idx_power_sat_applies_up | idx_power_sat_applies_dn;
        % eqns below from doi:10.1016/j.ifacol.2024.10.093
        alpha = imag(Z_th) ./ real(Z_th); % eqn 8
        epsilon = zeros(size(mag_X_f)); 
        epsilon(idx_force_sat_applies)  = 1; % effort limit
        epsilon((idx_power_sat_applies_up | idx_power_sat_applies_dn) & ~idx_force_sat_applies) = 1; % power limit phase is undefined, so treat it as an effort limit to minimize force for a given power limit.
        epsilon(idx_soln_applies & ~idx_force_sat_applies & ~(idx_power_sat_applies_up | idx_power_sat_applies_dn)) = -1; % flow limit
        sigma = sqrt( (alpha.^2 .* ctrl_mult_mag.^2 + 1).^2 + alpha.^2 .* (ctrl_mult_mag.^2 + 1).^2);
        acos_argument = -2 * alpha .* ctrl_mult_mag ./ sigma;
        acos_argument = max(-1, min(1, acos_argument));
        atan_denominator = sigma + epsilon .* alpha .* (1+ctrl_mult_mag).^2;
        small_denominator = abs(atan_denominator) < eps;
        atan_denominator(small_denominator) = eps .* sign(atan_denominator(small_denominator));
        atan_denominator(small_denominator & atan_denominator == 0) = eps;
        atan_argument = (alpha.^2 .*  ctrl_mult_mag.^2 + 1) ./ atan_denominator;
        ctrl_mult_phase_desired = 2 * atan(atan_argument) + epsilon .* acos(acos_argument); % eqn 4
        optimality_err(idx_soln_applies) = wrapToPi(ctrl_mult_phase(idx_soln_applies) - ctrl_mult_phase_desired(idx_soln_applies));

        P_unsat_safe = max(P_unsat, eps);
        power_lost = (P_unsat - P_sat) ./ P_unsat_safe;

        force_const_active = ismembertol(mag_U,   4/pi * F_max);
        min_f_const_active = ismembertol(mag_X_f, X_f_lower_limit);
        max_f_const_active = ismembertol(mag_X_f, X_f_upper_limit);
        min_s_const_active = ismembertol(mag_X_s, X_s_lower_limit);
        max_s_const_active = ismembertol(mag_X_s, X_s_upper_limit);
        max_u_const_active = ismembertol(mag_X_u, X_u_upper_limit);
        max_p_const_active = ismembertol(P_sat,   P_max);
        min_p_const_active = ismembertol(P_sat,   P_min);

        any_constr_active = force_const_active | min_f_const_active | max_f_const_active ...
                          | min_s_const_active | max_s_const_active | max_u_const_active ...
                          | max_p_const_active | min_p_const_active;
        any_constr_violated = constr_violation > 0 & ~any_constr_active;
        all_constr_inactive = ~any_constr_active & ~any_constr_violated;
        
        optimality_err(mult_const_viol_unsat & any_constr_active) = 0;
        optimality_err(mult_const_viol_unsat & any_constr_violated) = 10;
        optimality_err(mult_const_viol_unsat & all_constr_inactive) = ...
            power_lost(mult_const_viol_unsat & all_constr_inactive);
        
    end
end

function [X,angle_X] = second_order_transfer_fcn(w,m,b,k,F,F_phase)
    imag_term = b .* w;
    real_term = k - m .* w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    X = X_over_F_mag .* F;
    if nargout > 1
        angle_F_over_X = atan2(imag_term,real_term);
        angle_X = F_phase - angle_F_over_X;
    end
end

function ang = my_angle(complex_term)
    ang = atan2(imag(complex_term),real(complex_term));
end
