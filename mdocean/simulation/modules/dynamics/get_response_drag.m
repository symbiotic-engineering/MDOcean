function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         qcqp_debug,...
         F_drag_f,F_drag_s,...
    phase_F_drag_f,phase_F_drag_s] = get_response_drag(H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                        gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,...
                                        drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                        X_max,control_type,multibody,merge_bodies,...
                                        X_tol,phase_X_tol,F_lim_tol,X_lim_tol,...
                                        max_drag_iters_fixed_point,max_drag_iters_solver,...
                                        drag_convergence_plot_on,drag_fcn,...
                                        D_f,D_f_in,D_d,T_f_slam,T_s_slam)

    % initial guess: 2m float amplitude, 0.5m spar amplitude
    X_f_guess = 2 * ones(size(w));
    phase_X_f_guess = zeros(size(w));
    X_s_guess = .5 * ones(size(w));
    phase_X_s_guess = zeros(size(w));

    qcqp_debug = struct('centers', [], 'radii', [], 'labels', {{}}, ...
                        'Gamma_opt', NaN, 'alpha', NaN, 'Z_th', NaN, ...
                        'w', NaN, 'n_active_constraints', 0, 'feasible', false);

    % package inputs
    dynam_inputs = {H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,...
                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                    X_max,control_type,multibody,merge_bodies,...
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
         F_drag_f,F_drag_s,...
         phase_F_drag_f,phase_F_drag_s] = fixed_point_iteration(X_f_guess,X_s_guess,...
                                                    phase_X_f_guess,phase_X_s_guess,...
                                                    dynam_inputs{:}, ...
                                                    max_drag_iters_fixed_point);

        idx_use = isfinite(mag_X_f) | isnan(B_h_f);

        % update guesses
        [X_f_guess(idx_use),phase_X_f_guess(idx_use),...
         X_s_guess(idx_use),phase_X_s_guess(idx_use)] = deal(mag_X_f(idx_use),phase_X_f(idx_use),...
                                                       mag_X_s(idx_use),phase_X_s(idx_use));

    end
    % then do nonlinear solver to finish
    if max_drag_iters_solver > 0
        
        [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         qcqp_debug,...
         F_drag_f,F_drag_s,...
         phase_F_drag_f,phase_F_drag_s] = solver(X_f_guess,X_s_guess,...
                                        phase_X_f_guess,phase_X_s_guess,...
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
     qcqp_debug,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,...
     phase_F_drag_s] = solver(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                                F_max,drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody,merge_bodies,...
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

    % anonymous function to take 4 sea state matrices and turn into collapsed row vector
    flatten = @(X1,X2,X3,X4) [reshape(X1(idx_not_nan),[],1);...
                               reshape(X2(idx_not_nan),[],1);...
                               reshape(X3(idx_not_nan),[],1);...
                               reshape(X4(idx_not_nan),[],1)];

    % prepare inputs for solver
    x0 = flatten(X_f_guess, X_s_guess, phase_X_f_guess, phase_X_s_guess);

    % fun_inner takes 1 input and returns the 4+ outputs of dynamics_error_wrapper
    fun_inner = @(x) dynamics_error_wrapper(unflatten(x,1),unflatten(x,2),...
                                    unflatten(x,3),unflatten(x,4),...
                                    H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,...
                                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                    X_max,control_type,multibody,merge_bodies,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam);
    
    sparsity = repmat(eye(N_ss_nz),4);
    opts = optimoptions('fsolve','JacobPattern',sparsity,...
                                 'MaxIterations',max_drag_iters,...
                                 'Algorithm','trust-region',...
                                 ...%'StepTolerance',min([X_tol,phase_X_tol,F_lim_tol,X_lim_tol]),...
                                 'Display','off');
    if drag_convergence_plot_on
        opts.PlotFcn = {'optimplotx','optimplotfval'};
    end

    % solve
    [x_solved,~,~,~] = fsolve(@(x)fun_outer(x,fun_inner,flatten),x0,opts);

    % unpack
    [X_f_solved,X_s_solved,...
     phase_X_f_solved,phase_X_s_solved] = deal(unflatten(x_solved,1),...
                                               unflatten(x_solved,2),...
                                               unflatten(x_solved,3),...
                                               unflatten(x_solved,4));

    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     qcqp_debug,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,...
     phase_F_drag_s] = dynamics_from_guess(X_f_solved, phase_X_f_solved, mag_v0_f, drag_const_f, ...
                                    X_s_solved, phase_X_s_solved, mag_v0_s, drag_const_s, ...
                                    B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,H,w,k_wvn,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                    control_type,multibody,merge_bodies,F_max,X_max,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam);

end
function out_flat = fun_outer(x,fun_inner,flatten)
    [Y1,Y2,Y3,Y4] = fun_inner(x);
    out_flat = flatten(Y1,Y2,Y3,Y4);
end

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,...
     phase_F_drag_s] = fixed_point_iteration(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                                F_max,drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody,merge_bodies,...
                                                X_tol,phase_X_tol,F_lim_tol,X_lim_tol,...
                                                drag_convergence_plot_on,drag_fcn,...
                                                D_f, D_f_in, D_d, T_f_slam, ...
                                                T_s_slam, max_drag_iters)

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
         mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         F_drag_f,F_drag_s,...
         phase_F_drag_f,phase_F_drag_s] = dynamics_error_wrapper(X_f_guess,X_s_guess,...
                                                       phase_X_f_guess,phase_X_s_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,F_max,...
                                                drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody,merge_bodies,...
                                                drag_fcn, D_f, D_f_in, D_d, ...
                                                T_f_slam, T_s_slam);

        X_err       = max( abs([X_f_err,X_s_err]),             [], 'all');
        phase_X_err = max( abs([phase_X_f_err,phase_X_s_err]), [], 'all');

        % new guesses
        X_f_guess = mag_X_f;
        X_s_guess = mag_X_s;
        phase_X_f_guess = phase_X_f;
        phase_X_s_guess = phase_X_s;

        % check convergence (only drag convergence; QCQP is solved analytically)
        X_converged = X_err < X_tol;
        phase_X_converged = phase_X_err < phase_X_tol;
        converged = X_converged && phase_X_converged;

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
          mag_U,phase_U,...
          real_P,reactive_P,...
          mag_X_u,phase_X_u,...
          mag_X_f,phase_X_f,...
          mag_X_s,phase_X_s,...
          B_p,K_p,P_sat_ratio,...
          F_drag_f,F_drag_s,...
     phase_F_drag_f,phase_F_drag_s] = dynamics_error_wrapper(X_f_guess,X_s_guess,...
                                                        phase_X_f_guess,phase_X_s_guess,...
                                                H,w,k_wvn,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                                F_max,drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody,merge_bodies,drag_fcn,...
                                                D_f,D_f_in,D_d, T_f_slam, T_s_slam)
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     ~,...
     F_drag_f,F_drag_s,...
     phase_F_drag_f,phase_F_drag_s] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                    X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                    B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,H,w,k_wvn,...
                                    gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                    control_type,multibody,merge_bodies,F_max,X_max,...
                                    drag_fcn, D_f, D_f_in, D_d, T_f_slam, T_s_slam); % closed loop dynamics
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
    idx_unstable_or_fp_error = mag_X_f==Inf | (isnan(mag_X_f) & ~isnan(X_f_guess));
    X_f_err(idx_unstable_or_fp_error) = 0;
    X_s_err(idx_unstable_or_fp_error) = 0;
    phase_X_f_err(idx_unstable_or_fp_error) = 0;
    phase_X_s_err(idx_unstable_or_fp_error) = 0;
end

function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p_sat,K_p_sat,...
         P_sat_ratio,...
         qcqp_debug,...
     mag_F_drag_f,mag_F_drag_s,...
     phase_F_drag_f,phase_F_drag_s] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                        X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                        B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,H,w,k_wvn,...
                                        gamma_f_mag,gamma_f_phase,gamma_s_mag,gamma_s_phase,...
                                        control_type,multibody,merge_bodies,F_max,X_max,...
                                        drag_fcn,...
                                        D_f, D_f_in, D_d, T_f_slam, T_s_slam)

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



    if F_max==0
        [B_p,K_p] = deal(zeros(size(w)));
    else
        % unsaturated optimal control gains
        [B_p,K_p] = controller(real_G_u, imag_G_u, w, control_type);
    end
    % unsaturated response (stabilized) - capture phases for QCQP
    stabilize_B = true;
    stabilize_K = strcmpi(control_type,'reactive');
    [mag_U_unsat,phase_U_unsat,P_unsat,~,...
     mag_X_u_unsat,phase_X_u_unsat,...
     mag_X_f_unsat,phase_X_f_unsat,...
     mag_X_s_unsat,phase_X_s_unsat,...
     B_p_stabilized,K_p_stabilized] = control_evaluation_fcn(K_p,B_p,stabilize_B,stabilize_K);
    
    Z_th = 1 ./ (real_G_u + 1i * imag_G_u);

    % solve constrained optimal control via QCQP (circle intersection)
    [B_p_sat,K_p_sat,qcqp_debug] = solve_qcqp_control(Z_th, w, ...
                            mag_X_u_unsat, phase_X_u_unsat, ...
                            mag_X_f_unsat, phase_X_f_unsat, ...
                            mag_X_s_unsat, phase_X_s_unsat, ...
                            mag_U_unsat, phase_U_unsat, ...
                            B_p_stabilized, K_p_stabilized, ...
                            F_max, X_max, ...
                            B_h_f, B_h_s, B_c, K_h_f, K_h_s, ...
                            B_drag_f, B_drag_s, m_c, m_f, m_s, ...
                            control_type, multibody, merge_bodies);

    % evaluate final response with constrained optimal controller
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);

    P_sat_ratio = real_P ./ P_unsat;
    F_drag_f = -B_drag_f .* mag_X_f .* w .* exp(1i * (phase_X_f + pi/2)) + gamma_drag_f .* H/2;
    F_drag_s = -B_drag_s .* mag_X_s .* w .* exp(1i * (phase_X_s + pi/2)) + gamma_drag_s .* H/2;
    mag_F_drag_f = abs(F_drag_f);
    mag_F_drag_s = abs(F_drag_s);
    phase_F_drag_f = angle(F_drag_f);
    phase_F_drag_s = angle(F_drag_s);
end

function [B_p_sat,K_p_sat,qcqp_debug] = solve_qcqp_control(Z_th, w, ...
                            mag_X_u_unsat, phase_X_u_unsat, ...
                            mag_X_f_unsat, phase_X_f_unsat, ...
                            mag_X_s_unsat, phase_X_s_unsat, ...
                            mag_U_unsat, phase_U_unsat, ...
                            B_p_stabilized, K_p_stabilized, ...
                            F_max, X_max, ...
                            B_h_f, B_h_s, B_c, K_h_f, K_h_s, ...
                            B_drag_f, B_drag_s, m_c, m_f, m_s, ...
                            control_type, multibody, merge_bodies)
% Solve the constrained optimal control problem analytically using
% the QCQP circle intersection method from the paper.
% Returns the constrained optimal B_p_sat and K_p_sat.

    X_f_max = X_max(1);
    X_s_max = X_max(2);
    X_u_max = X_max(3);
    
    % small tolerance for negligible coefficients (dimensionless and dimensional)
    COEFF_TOL = 1e-10;
    GAMMA_TOL = 1e-12;
    
    B_p_sat = B_p_stabilized;
    K_p_sat = K_p_stabilized;
    qcqp_debug = struct('centers', [], 'radii', [], 'labels', {{}}, ...
                        'Gamma_opt', NaN, 'alpha', NaN, 'Z_th', NaN, ...
                        'w', NaN, 'n_active_constraints', 0, 'feasible', false);
    
    % skip QCQP if no constraints or no PTO
    use_force_sat = isfinite(F_max) && F_max > 0;
    use_amp_sat = any(isfinite(X_max));
    if ~use_force_sat && ~use_amp_sat
        return
    end
    if F_max == 0
        return
    end
    
    % compute complex phasors for the unsaturated (Gamma=0) response
    % and the derivative coefficients c_Y for each quantity Y
    % Y(Gamma) = Y_0 + c_Y * Gamma, where Y_0 is the unsaturated value
    
    sz = size(w);
    
    for idx = 1:numel(w)
        if isnan(B_h_f(idx))
            continue  % skip NaN sea states
        end
        
        Z_th_i = Z_th(idx);
        w_i = w(idx);
        
        if ~isfinite(Z_th_i) || real(Z_th_i) <= 0
            continue  % skip unstabilizable sea states
        end
        
        % complex phasors at Gamma=0
        X_u_0 = mag_X_u_unsat(idx) * exp(1i * phase_X_u_unsat(idx));
        X_f_0 = mag_X_f_unsat(idx) * exp(1i * phase_X_f_unsat(idx));
        X_s_0 = mag_X_s_unsat(idx) * exp(1i * phase_X_s_unsat(idx));
        U_0   = mag_U_unsat(idx)   * exp(1i * phase_U_unsat(idx));
        
        % I_p = complex current phasor at Gamma=0
        % I = I_p*(1-Gamma), so at Gamma=0: I_0 = I_p
        % velocity = i*w*X_u, so I_p = i*w_i*X_u_0
        I_p = 1i * w_i * X_u_0;
        
        if abs(I_p) < COEFF_TOL
            continue  % no response, nothing to constrain
        end
        
        % derivative coefficients: Y(Gamma) = Y_0 + c_Y * Gamma
        % From the paper: V = I_p*Z_th^**(1+Gamma), I = I_p*(1-Gamma)
        c_V = I_p * conj(Z_th_i);   % force derivative: dV/dGamma
        c_I = -I_p;                   % velocity derivative: dI/dGamma
        c_Xu = -X_u_0;               % PTO displacement: X_u = I/(iw)
        
        % float and spar displacement derivatives depend on topology
        if ~multibody || merge_bodies
            % single body: X_f = X_u, X_s = 0
            c_Xf = c_Xu;
            c_Xs = 0;
        else
            % multibody: compute transfer functions from impedance matrices
            % All arrays indexed by (idx) to extract the scalar for this sea state
            B_f_total = B_h_f(idx) + B_drag_f(idx);
            B_s_total = B_h_s(idx) + B_drag_s(idx);
            Z_f_i = B_f_total + 1i*(m_f(idx)*w_i - K_h_f(idx)/w_i);
            Z_s_i = B_s_total + 1i*(m_s(idx)*w_i - K_h_s(idx)/w_i);
            Z_c_i = B_c(idx) + 1i*(m_c(idx)*w_i);
            det_Z_i = Z_f_i .* Z_s_i - Z_c_i.^2;
            
            % transfer from PTO force to body displacement
            % v_f = v_f_forced + (Z_s+Z_c)/det_Z * F_pto
            % X_f = X_f_forced + (Z_s+Z_c)/(det_Z*iw) * F_pto
            % F_pto = V = I_p*Z_th^**(1+Gamma)
            % c_Xf = (Z_s+Z_c)/(det_Z*iw) * I_p * Z_th^*
            T_f = (Z_s_i + Z_c_i) ./ (det_Z_i .* 1i .* w_i);
            T_s = -(Z_f_i + Z_c_i) ./ (det_Z_i .* 1i .* w_i);
            c_Xf = T_f .* I_p .* conj(Z_th_i);
            c_Xs = T_s .* I_p .* conj(Z_th_i);
        end
        
        % build circle constraints
        % For constraint |Y| <= Y_max where Y = Y_0 + c_Y * Gamma:
        % circle center = -Y_0/c_Y, radius = Y_max/|c_Y|, S=+1 (inside)
        centers = [];
        radii = [];
        current_labels = {};
        
        % force constraint: |V| <= F_max
        if use_force_sat && abs(c_V) > COEFF_TOL
            center_F = -U_0 / c_V;
            radius_F = F_max / abs(c_V);
            if mag_U_unsat(idx) > F_max  % only add if violated
                centers = [centers; real(center_F), imag(center_F)];
                radii = [radii; radius_F];
                current_labels{end+1} = 'Force limit';
            end
        end
        
        % float amplitude constraint: |X_f| <= X_f_max
        if isfinite(X_f_max) && abs(c_Xf) > COEFF_TOL
            center_Xf = -X_f_0 / c_Xf;
            radius_Xf = X_f_max / abs(c_Xf);
            if mag_X_f_unsat(idx) > X_f_max  % only add if violated
                centers = [centers; real(center_Xf), imag(center_Xf)];
                radii = [radii; radius_Xf];
                current_labels{end+1} = 'Float amplitude';
            end
        end
        
        % spar amplitude constraint: |X_s| <= X_s_max
        if isfinite(X_s_max) && abs(c_Xs) > COEFF_TOL && multibody && ~merge_bodies
            center_Xs = -X_s_0 / c_Xs;
            radius_Xs = X_s_max / abs(c_Xs);
            if mag_X_s_unsat(idx) > X_s_max  % only add if violated
                centers = [centers; real(center_Xs), imag(center_Xs)];
                radii = [radii; radius_Xs];
                current_labels{end+1} = 'Spar amplitude';
            end
        end
        
        % PTO amplitude constraint: |X_u| <= X_u_max
        if isfinite(X_u_max) && abs(c_Xu) > COEFF_TOL
            center_Xu = -X_u_0 / c_Xu;
            radius_Xu = X_u_max / abs(c_Xu);
            if mag_X_u_unsat(idx) > X_u_max  % only add if violated
                centers = [centers; real(center_Xu), imag(center_Xu)];
                radii = [radii; radius_Xu];
                current_labels{end+1} = 'PTO amplitude';
            end
        end
        
        % positive power constraint: |Gamma - i*alpha|^2 <= 1 + alpha^2
        % where alpha = Im(Z_th)/Re(Z_th)
        alpha = imag(Z_th_i) / real(Z_th_i);
        center_P = [0, alpha];
        radius_P = sqrt(1 + alpha^2);
        % always include positive power constraint (it's cheap and prevents P<0)
        centers = [centers; center_P];
        radii = [radii; radius_P];
        current_labels{end+1} = 'Positive power';
        
        % damping control: add Q=0 constraint circle
        % Gamma must lie ON the Q=0 circle (equality constraint)
        % Q=0 circle: center = -i*R/X, radius = |Z_th|/|X|
        % where R = Re(Z_th), X = Im(Z_th)
        p_star = [];
        if strcmpi(control_type, 'damping') && abs(imag(Z_th_i)) > COEFF_TOL
            R_th = real(Z_th_i);
            X_th = imag(Z_th_i);
            center_Q0 = [0, -R_th/X_th];
            radius_Q0 = abs(Z_th_i) / abs(X_th);
            % for damping, find optimal on Q=0 circle subject to other constraints
            Gamma_opt = solve_damping_qcqp(center_Q0, radius_Q0, centers, radii);
        else
            % reactive control or no reactance: use standard circle intersection
            if isempty(centers)
                Gamma_opt = 0;  % no constraints active
            else
                [p_star, ~, ~] = circle_intersect_optim(centers, radii);
                if isempty(p_star)
                    % Constraint circles are mutually infeasible.
                    % This typically happens when force (circle at Gamma=-1)
                    % and amplitude (circle at Gamma=+1) limits are both
                    % violated and r_F + r_Xu < 2 (they can't intersect).
                    % Fall back: try each individual constraint circle with
                    % the positive-power circle to find best feasible Gamma.
                    % Power circle is always the last entry in centers/radii.
                    Gamma_opt = 1;  % default: shutdown
                    best_norm = Inf;
                    n_c = size(centers, 1);
                    c_pow = centers(n_c, :);
                    r_pow = radii(n_c);
                    for k = 1:n_c - 1
                        % closest boundary point on circle k to origin
                        dist_to_center_k = norm(centers(k,:));
                        if dist_to_center_k < COEFF_TOL
                            boundary_point_k = [radii(k), 0];
                        else
                            boundary_point_k = centers(k,:) - radii(k) * centers(k,:)/dist_to_center_k;
                        end
                        % accept only if inside the positive-power circle
                        if norm(boundary_point_k - c_pow) <= r_pow + 1e-4 && norm(boundary_point_k) < best_norm
                            best_norm = norm(boundary_point_k);
                            Gamma_opt = boundary_point_k(1) + 1i*boundary_point_k(2);
                        end
                    end
                else
                    Gamma_opt = p_star(1) + 1i*p_star(2);
                end
            end
        end
        
        % update debug struct if this is the most constrained feasible sea state
        n_non_power = size(centers, 1) - 1;  % last circle is always the power circle
        if n_non_power > qcqp_debug.n_active_constraints && ~isempty(centers) && isfinite(Gamma_opt)
            qcqp_debug.centers = centers;
            qcqp_debug.radii = radii;
            qcqp_debug.labels = current_labels;
            qcqp_debug.Gamma_opt = Gamma_opt;
            qcqp_debug.alpha = imag(Z_th_i) / real(Z_th_i);
            qcqp_debug.Z_th = Z_th_i;
            qcqp_debug.w = w_i;
            qcqp_debug.n_active_constraints = n_non_power;
            qcqp_debug.feasible = ~isempty(p_star);
        end

        % convert Gamma to Z_l, then to B_p and K_p
        if abs(Gamma_opt) < GAMMA_TOL
            % no constraint active, keep unsaturated gains
            continue
        end
        
        Z_l = conj(Z_th_i) * (1 + Gamma_opt) / (1 - Gamma_opt);
        B_p_sat(idx) = real(Z_l);
        K_p_sat(idx) = -w_i * imag(Z_l);
        
        % ensure B_p >= 0 (stabilize)
        if B_p_sat(idx) < 0
            B_p_sat(idx) = 0;
        end
    end
end

function Gamma_opt = solve_damping_qcqp(center_Q0, radius_Q0, centers, radii)
% For damping control, find the point on the Q=0 circle closest to the 
% origin that satisfies all other constraint circles.
% The Q=0 circle parametrically: Gamma(t) = center + radius * exp(i*t)

    % closest point on Q=0 circle to origin (unconstrained damping optimal)
    c_Q0 = center_Q0(1) + 1i*center_Q0(2);
    dist_to_center = abs(c_Q0);
    if dist_to_center < 1e-10
        % center at origin, any point on circle is equally close
        Gamma_closest = radius_Q0;
    else
        Gamma_closest = c_Q0 - radius_Q0 * c_Q0 / dist_to_center;
    end
    
    % check if the closest point is feasible
    p_closest = [real(Gamma_closest), imag(Gamma_closest)];
    if isempty(centers) || check_inside_circles(p_closest, centers, radii)
        Gamma_opt = Gamma_closest;
        return
    end
    
    % find intersections of Q=0 circle with each constraint circle
    candidates = [];
    for i = 1:length(radii)
        [xs, ys] = circle_circle_intersect_inline(...
            center_Q0(1), center_Q0(2), radius_Q0, ...
            centers(i,1), centers(i,2), radii(i));
        if ~isnan(xs(1))
            candidates = [candidates; xs(1), ys(1); xs(2), ys(2)]; %#ok<AGROW>
        end
    end
    
    % filter feasible candidates
    feasible = [];
    for i = 1:size(candidates,1)
        p = candidates(i,:);
        % check on Q=0 circle and inside all other constraint circles
        on_Q0 = abs(norm(p - center_Q0) - radius_Q0) < 1e-4;
        if on_Q0 && check_inside_circles(p, centers, radii)
            feasible = [feasible; p]; %#ok<AGROW>
        end
    end
    
    if isempty(feasible)
        % No point on the Q=0 circle satisfies all constraints simultaneously.
        % Fall back to the unconstrained damping-optimal point (Gamma_closest),
        % which always satisfies the positive-power constraint (the Q=0 and
        % power circles always intersect when Im(Z_th)≠0).
        % This may violate force/amplitude limits but avoids complete shutdown.
        Gamma_opt = Gamma_closest;
    else
        % choose closest to origin
        dists = vecnorm(feasible, 2, 2);
        [~, best] = min(dists);
        Gamma_opt = feasible(best,1) + 1i*feasible(best,2);
    end
end

function inside = check_inside_circles(p, centers, radii)
% Check if point p=[x,y] is inside all circles
    dists = vecnorm(p - centers, 2, 2);
    inside = all(dists <= radii + 1e-4);
end

function [xs, ys] = circle_circle_intersect_inline(x1,y1,r1,x2,y2,r2)
% Inline circle-circle intersection (same as in circle_intersect_optim)
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
    if d > r1+r2 || d < abs(r1-r2) || d == 0
        xs = [NaN, NaN];
        ys = [NaN, NaN];
        return
    end
    a = (r1^2 - r2^2 + d^2) / (2*d);
    h = sqrt(max(r1^2 - a^2, 0));
    mx = x1 + a*(x2-x1)/d;
    my = y1 + a*(y2-y1)/d;
    dx = h*(y2-y1)/d;
    dy = h*(x2-x1)/d;
    xs = [mx+dx, mx-dx];
    ys = [my-dy, my+dy];
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
    F_f_phase = angle(F_f);
    F_s_phase = angle(F_s);

    % thevenin admittance of the controlled DOF as seen from load
    if multibody
        [real_G_u,imag_G_u] = multibody_impedance(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w);
        
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
    end

    idx_not_stabilizable = real_G_u < 0; % controlled dof sees negative source 
                                         % impedance, so it's not only unstable
                                         % but also not stabilizable


    % function to evaluate any controller for fixed linear system
    control_evaluation_fcn = @(K_p,B_p,stabilize_B,stabilize_K) evaluate_any_controller(...
                                        B_c, B_f, B_s, K_f, K_s, ...
                                        m_c, m_f, m_s, w, K_p, B_p, ...
                                        F_f_mag, F_f_phase, F_s_mag, ...
                                        F_s_phase, multibody, merge_bodies,...
                                        idx_not_stabilizable, ...
                                        stabilize_B, stabilize_K);
end

function [B_p,K_p] = controller(real_G_u, imag_G_u, w, control_type)

    mag_G_u_squared = real_G_u.^2 + imag_G_u.^2;

    % control - set powertrain coefficients
    if strcmpi(control_type,'reactive')
        K_p = -w .* imag_G_u ./ mag_G_u_squared;
        B_p = real_G_u ./ mag_G_u_squared;
    elseif strcmpi(control_type, 'damping')
        B_p = 1 ./ sqrt(mag_G_u_squared);
        K_p = 1e-8; % can't be quite zero because r_k = Inf
    end

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
    G_integral_weighted = G_integral_real_weighted + 1i*G_integral_imag_weighted;

    % excitation only when mag_v = 0
    mag_v_or_v0 = mag_v;
    mag_v_or_v0(mag_v == 0) = mag_v0(mag_v == 0);

    % derived p141 of notebook 9 3/5/26
    mag_v_constant_term = drag_const / (pi*(1-alpha^2)) * 2 * mag_v_or_v0;
    B_drag_2   = mag_v_constant_term           .* B_integral_weighted;
    assert(all( B_drag_2(isfinite(B_drag_2)) >= 0 ),'negative damping from integral method')
    gamma_drag = mag_v_constant_term .* (mag_v0 ./ (H/2)) .* G_integral_weighted;

    plot_on = false;
    drag_debug = false; % set true to plot new integral drag vs old drag
    if drag_debug
        % uncomment the following to plot once solver has converged
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

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_l,K_l] = evaluate_any_controller(...
                            B_c, B_f, B_s, K_f, K_s, ...
                            m_c, m_f, m_s, w, K_l, B_l, ...
                            F_f_mag, F_f_phase, F_s_mag, F_s_phase, ...
                            multibody, merge_bodies, idx_not_stabilizable,...
                            stabilize_B,stabilize_K)

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
        F_m_phase = angle(F_exc);

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
        phase_U_f = angle(F_brake_f);
        mag_U_s = abs(F_brake_s);
        phase_U_s = angle(F_brake_s);
        % fixme it should be mag_U = mag_U_s and phase_U = pi+phase_U_s,
        % take the max and use float phase for now
        mag_U = max(mag_U_f,mag_U_s);
        phase_U = phase_U_f;

    else
        if multibody
            [mag_U,phase_U,...
                real_P,reactive_P,...
                mag_X_u,phase_X_u,...
                mag_X_f,phase_X_f,...
                mag_X_s,phase_X_s] = multibody_response(B_c, B_f, B_s, K_f, K_s, ...
                                                        m_c, m_f, m_s, w, ...
                                                        K_l, B_l, ...
                                                        F_f_mag, F_f_phase, ...
                                                        F_s_mag, F_s_phase);
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
            reactive_P = 0; % fixme this is incorrect but doesn't affect anything rn
        end
    end

    % set values to Inf when closed loop unstable
    mag_U(idx_closed_loop_unstable)      = Inf;
    phase_U(idx_closed_loop_unstable)    = Inf;
    real_P(idx_closed_loop_unstable)     = Inf;
    reactive_P(idx_closed_loop_unstable) = Inf;
    mag_X_u(idx_closed_loop_unstable)    = Inf;
    phase_X_u(idx_closed_loop_unstable)  = Inf;
    mag_X_f(idx_closed_loop_unstable)    = Inf;
    phase_X_f(idx_closed_loop_unstable)  = Inf;
    mag_X_s(idx_closed_loop_unstable)    = Inf;
    phase_X_s(idx_closed_loop_unstable)  = Inf;
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

    Z_l = B_l_sat - 1i*K_l_sat./w;
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

function [X,angle_X] = second_order_transfer_fcn(w,m,b,k,F,F_phase)
    imag_term = b .* w;
    real_term = k - m .* w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    X = X_over_F_mag .* F;
    if nargout > 1
        X_over_F_phase = atan2(imag_term,real_term);
        angle_X = X_over_F_phase + F_phase;
    end
end

