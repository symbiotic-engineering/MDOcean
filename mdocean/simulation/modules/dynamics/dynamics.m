function [F_heave_storm, F_surge_storm, ...
    F_heave_op, F_surge_op, F_ptrain_max, ...
    P_var, P_avg_elec, P_matrix_elec, ...
    X_constraints, B_p, K_p, mag_U, ...
    X_u, X_f, X_s, ...
    P_matrix_mech,P_sat_ratio,...
    A_f_over_rho, A_s_over_rho, A_c_over_rho,...
    B_f_over_rho_w, B_s_over_rho_w,B_c_over_rho_w,...
    gamma_f_over_rho_g, gamma_s_over_rho_g,...
    gamma_phase_f,gamma_phase_s,w,phase_X_f,phase_X_u] = dynamics(in,m_float,m_spar)

    % use probabilistic sea states for power and PTO force and max amplitude
    [T,Hs] = meshgrid(in.T,in.Hs);
    zero_prob_idxs = in.JPD==0;
    T(zero_prob_idxs) = NaN;
    Hs(zero_prob_idxs) = NaN;

    if in.use_force_sat
        F_limit = in.F_max;
    else
        F_limit = Inf;
    end

    T_f_slam_op    = min(in.T_f_1,in.h_f-in.T_f_2);
    T_f_slam_storm = min(in.T_f_2,in.h_f-in.T_f_2);

    [P_matrix_mech,X_constraints_op,B_p,...
    K_p,mag_U,X_u,X_f,X_s,F_heave_op,...
    F_surge_op,F_ptrain_max,P_sat_ratio,...
    A_f_over_rho,A_s_over_rho, A_c_over_rho,...
    B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
    gamma_f_over_rho_g, gamma_s_over_rho_g, ...
    gamma_phase_f, gamma_phase_s,w,...
    phase_X_f,phase_X_u] = get_power_force(in,T,Hs,m_float,m_spar,...
                                            F_limit, zero_prob_idxs, ...
                                            T_f_slam_op, in.power_scale_coeffs,...
                                            in.operational_drag_convergence_plot_on);
    
    % account for powertrain electrical losses
    P_matrix_elec = P_matrix_mech * in.eff_pto;
    
    % saturate maximum power
    if in.use_power_sat
        P_matrix_elec = min(P_matrix_elec,in.P_max,'includenan');
    end
    
    % weight power across all sea states
    P_weighted = P_matrix_elec .* in.JPD / 100;
    P_avg_elec = sum(P_weighted(:),'omitnan'); 
    
    assert(isreal(P_avg_elec))
    
    % use max sea states for structural forces
    in.use_multibody = false; % storm is alway singlebody dynamics because brake clamps float to spar
    [~,X_constraints_storm,~,~,~,~,~,~,F_heave_storm,F_surge_storm] = get_power_force(in, ...
                                in.T_struct, in.Hs_struct, m_float, m_spar, ...
                                0, false(size(in.T_struct)), T_f_slam_storm, [1 0 0 1], false);
    F_heave_storm = F_heave_storm * in.F_heave_mult;
    
    % use all X constraints operationally, only use slamming in storm
    X_constraints_storm = X_constraints_storm(6:end);
    X_constraints_storm = 1 + 0*X_constraints_storm; % fixme this overrides the constraint
    X_constraints = [X_constraints_op X_constraints_storm];

    % coefficient of variance (normalized standard deviation) of power
    P_var = std(P_matrix_elec(:), in.JPD(:),'omitnan') / P_avg_elec;
    P_var = P_var * 100; % convert to percentage
    

end

function [P_matrix, X_constraints, B_p, K_p, mag_U, mag_X_u, mag_X_f, mag_X_s,...
          F_heave_f_max, F_surge, F_ptrain_max, P_sat_ratio,...
          A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
        B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
        gamma_f_over_rho_g, gamma_s_over_rho_g, ...
        gamma_phase_f, gamma_phase_s,w,phase_X_f,phase_X_u] = get_power_force(in,T,Hs, m_float, m_spar, ...
                                                            F_max, idx_constraint, T_f_slam, power_scale_coeffs,...
                                                            drag_convergence_plot_on)

    % get dynamic coefficients for float and spar
    % fixme: eventually should use in.D_f_in to allow a radial gap between float and spar
    [m_f,B_h_f,K_h_f,F_f_mag,F_f_phase,...
     m_s,B_h_s,K_h_s,F_s_mag,F_s_phase,...
     m_c,B_c,drag_const_f,drag_const_s,...
     mag_v0_f,mag_v0_s,w,k_wvn,A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
        B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
        gamma_f_over_rho_g, gamma_s_over_rho_g, ...
        gamma_phase_f, gamma_phase_s] = get_dynamic_coeffs(Hs, T, ...
                                            in.D_f, in.T_f_2, in.D_s, in.D_d, in.T_s, in.h, ...
                                            m_float, m_spar, in.spar_excitation_coeffs,...
                                            in.C_d_float, in.C_d_spar, ...
                                            in.rho_w, in.g, ...
                                            in.use_MEEM, in.harmonics, in.hydro);

    X_u_max = 1e6;%min(Hs / (2*sqrt(2)), in.T_f);

    % get response: includes drag and force saturation
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio] = get_response_drag(w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                            F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                            drag_const_f,drag_const_s,mag_v0_f,mag_v0_s, ...
                                            X_u_max,in.control_type,in.use_multibody,...
                                            in.X_tol,in.phase_X_tol,in.max_drag_iters,...
                                            drag_convergence_plot_on);

% FIXME: check stability of closed loop multibody system
    
    % apply empirical fitted freq-dependent scale factor to adjust power from 
    % singlebody to multibody, to match RM3 report
    if in.use_multibody==false
        c = power_scale_coeffs;
        power_scale_factor = c(1) ./ ( c(2) * T.^2 + c(3) * T + c(4) );
        P_matrix = real_P .* power_scale_factor;
    else
        P_matrix = real_P * in.power_scale_multibody;
    end
    
    % set values where JPD=0 to 0 to not include them in constraint
    mag_X_u_const = mag_X_u;
    mag_X_u_const(idx_constraint) = 0;
    mag_X_f_const = mag_X_f;
    mag_X_f_const(idx_constraint) = 0;
    mag_X_s_const = mag_X_s;
    mag_X_s_const(idx_constraint) = 0;

    X_u_max = max(mag_X_u_const,[],'all');
    X_f_max = max(mag_X_f_const,[],'all');
    X_s_max = max(mag_X_s_const,[],'all');

    % extra height on spar after accommodating relative displacement
    h_s_extra_up = (in.h_s - in.T_s - (in.h_f - in.T_f_2) - X_u_max) / in.h_s;
    h_s_extra_down = (in.T_s - in.T_f_2 - in.h_d - X_u_max) / in.h_s;

    % sufficient length of float support tube
    h_fs_extra = in.h_fs_clear / X_u_max - 1;

    % prevent violation of linear wave theory - notebook p153 9/20/25
    X_max_linear_f = 1/10 * (in.h - in.T_f_2);
    X_max_linear_s = 1/10 * (in.h - in.T_s);
    
    X_below_linear_f = X_max_linear_f / X_f_max - 1;
    X_below_linear_s = X_max_linear_s / X_s_max - 1;

    % prevent rising out of the water (slamming)
    wave_amp = Hs/(2*sqrt(2));
    theta_slam = max(0, -k_wvn * in.D_f / 2 + abs(pi - phase_X_f));
    X_slam = sqrt( T_f_slam^2 - (wave_amp .* sin(theta_slam)).^2 ) - wave_amp .* cos(theta_slam);
    X_slam( imag(X_slam)~=0 ) = 0; % case where slamming occurs even for stationary body
    X_below_wave = X_slam ./ mag_X_f_const - 1;

    plot_slamming = false;
    if plot_slamming
        X_star = (X_slam - T_f_slam)./wave_amp;
        X_slam_simple = T_f_slam - wave_amp;
        X_below_simple = X_slam_simple ./ mag_X_u_const - 1;
        make_slamming_plot(T,Hs,theta_slam,X_below_wave,X_star,X_below_simple)
    end

    X_below_wave(~isfinite(X_below_wave)) = 1; % constraint always satisfied when JPD=0

    X_constraints = [h_s_extra_up, h_s_extra_down, h_fs_extra, X_below_linear_f, X_below_linear_s, X_below_wave(:).'];

    % calculate forces
    if nargout > 3
        % set values where JPD=0 to 0 to not include them in constraint
        mag_U_const = mag_U;
        mag_U_const(idx_constraint) = 0;
            
        % powertrain force
        F_ptrain_max = max(mag_U_const,[],'all');
        F_ptrain_max = min(F_ptrain_max, F_max);

        % heave force: includes powertrain force and D'Alembert force        
        F_heave_f = combine_ptrain_dalembert_forces(m_float, w, mag_X_f_const, phase_X_f, mag_U_const, phase_U, F_max);
        F_heave_s = combine_ptrain_dalembert_forces(m_spar,  w, mag_X_s_const, phase_X_s, mag_U_const, phase_U, F_max);

        F_heave_f_max = max(F_heave_f,[],'all');
        F_heave_s_max = max(F_heave_s,[],'all');

        % surge force - from Eq 25 Newman 1963 - assumes slender kR << 1
        % https://apps.dtic.mil/sti/tr/pdf/AD0406333.pdf  
        F_surge_coeff = 2 * in.rho_w * w.^2 .* wave_amp ./ k_wvn;
        F_surge_f = F_surge_coeff * pi * in.D_f^2/4 .* (       1             - exp(-k_wvn*in.T_f_2));
        F_surge_s = F_surge_coeff * pi * in.D_s^2/4 .* (exp(-k_wvn*in.T_f_2) - exp(-k_wvn*in.T_s));

        F_surge_f_max = max(F_surge_f(~idx_constraint),[],'all');
        F_surge_s_max = max(F_surge_s(~idx_constraint),[],'all');

        F_surge = [F_surge_f_max F_surge_s_max 0];

        if ~isfinite(F_heave_f_max)
            disp('F heave f:')
            disp(F_heave_f)
            disp('m float: ')
            disp(m_float)
            disp('w: ')
            disp(w)
            disp('mag X f const: ')
            disp(mag_X_f_const)
            disp('phase X f: ')
            disp(phase_X_f)
            disp('mag U const:')
            disp(mag_U_const)
            disp('phase U: ')
            disp(phase_U)
            disp('F max: ')
            disp(F_max)
        end

    end
end

function F_heave = combine_ptrain_dalembert_forces(mass, w, mag_X, phase_X, mag_U, phase_U, F_max)
    F_dAlembert_mag = mass .* w.^2 .* mag_X;
    F_dAlembert_phase = phase_X;

    F_heave_real = F_dAlembert_mag .* cos(F_dAlembert_phase) + mag_U .* cos(phase_U);
    F_heave_imag = F_dAlembert_mag .* sin(F_dAlembert_phase) + mag_U .* sin(phase_U);

    F_heave_fund = sqrt(F_heave_real.^2 + F_heave_imag.^2);

    F_heave = min(F_heave_fund, F_max + F_dAlembert_mag);

end

function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio] = get_response_drag(w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                        F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                        drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                        X_max,control_type,multibody,...
                                        X_tol,phase_X_tol,max_drag_iters,...
                                        drag_convergence_plot_on)
    % initial guess: 2m float amplitude, 0.5m spar amplitude, unsaturated
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

    max_drag_iters_fixed_point = 3;
    max_drag_iters_solver      = 200;%max_drag_iters;

    % package inputs
    dynam_inputs = {w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                    F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                    X_max,control_type,multibody,...
                    X_tol,phase_X_tol,drag_convergence_plot_on};

    % first do fixed point iteration
    if max_drag_iters_fixed_point > 0
        [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio,...
         ctrl_mult,phase_ctrl_mult] = fixed_point_iteration(X_f_guess,X_s_guess,...
                                                    phase_X_f_guess,phase_X_s_guess,...
                                                    ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                    dynam_inputs{:}, ...
                                                    max_drag_iters_fixed_point);

        idx_use = (ctrl_mult >= 0 & ctrl_mult < 1e4) | isnan(ctrl_mult);

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
         B_p,K_p,P_sat_ratio] = solver(X_f_guess,X_s_guess,...
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
     B_p,K_p,P_sat_ratio] = solver(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                                    ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                                drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody,...
                                                X_tol,phase_X_tol,drag_convergence_plot_on,max_drag_iters)

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

    % anonymous function to take 6 sea state matrices and turn into collapsed row vector
    flatten = @(X1,X2,X3,X4,X5,X6) [reshape(X1(idx_not_nan),[],1);...
                                    reshape(X2(idx_not_nan),[],1);...
                                    reshape(X3(idx_not_nan),[],1);...
                                    reshape(X4(idx_not_nan),[],1);...
                                    reshape(X5(idx_not_nan),[],1);...
                                    reshape(X6(idx_not_nan),[],1)    ];

    % prepare inputs for solver
    x0 = flatten(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                ctrl_mult_guess,phase_ctrl_mult_guess);

    % fun_inner takes 1 input and returns the 6+ outputs of dynamics_error_wrapper
    fun_inner = @(x) dynamics_error_wrapper(unflatten(x,1),unflatten(x,2),...
                                    unflatten(x,3),unflatten(x,4),...
                                    unflatten(x,5),unflatten(x,6),...
                                    w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                    F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                    drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                    X_max,control_type,multibody);
    
    sparsity = repmat(eye(N_ss_nz),6);
    opts = optimoptions('fsolve','JacobPattern',sparsity,...
                                 'MaxIterations',max_drag_iters,...
                                 'Algorithm','trust-region',...
                                 ...%'StepTolerance',min(X_tol,phase_X_tol),...
                                 'Display','off');
    if drag_convergence_plot_on
        opts.PlotFcn = {'optimplotx','optimplotfval'};
    end

    % solve
    [x_solved,err,flag,out] = fsolve(@(x)fun_outer(x,fun_inner,flatten),x0,opts);

    % unpack
    [X_f_solved,X_s_solved,...
     phase_X_f_solved,phase_X_s_solved,...
     ctrl_mult_solved,phase_ctrl_mult_solved] = deal(unflatten(x_solved,1),...
                                                     unflatten(x_solved,2),...
                                                     unflatten(x_solved,3),...
                                                     unflatten(x_solved,4),...
                                                     unflatten(x_solved,5),...
                                                     unflatten(x_solved,6));

    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio] = dynamics_from_guess(X_f_solved, phase_X_f_solved, mag_v0_f, drag_const_f, ...
                                    X_s_solved, phase_X_s_solved, mag_v0_s, drag_const_s, ...
                                    B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,w,...
                                    F_f_mag,F_f_phase,F_s_mag,F_s_phase,...
                                    control_type,multibody,F_max,X_max,...
                                    ctrl_mult_solved,phase_ctrl_mult_solved);

end
function out_flat = fun_outer(x,fun_inner,flatten)
    [Y1,Y2,Y3,Y4,Y5,Y6] = fun_inner(x);
    out_flat = flatten(Y1,Y2,Y3,Y4,Y5,Y6);
end

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     ctrl_mult_guess,...
     phase_ctrl_mult_guess] = fixed_point_iteration(X_f_guess,X_s_guess,phase_X_f_guess,phase_X_s_guess,...
                                                ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                                drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody,...
                                                X_tol,phase_X_tol,drag_convergence_plot_on,max_drag_iters)

    converged = false;
    iters = 0;

    % loop until converged
    while ~converged
        if drag_convergence_plot_on
            % record guesses
            % 10,1 was found manually to be the problem sea state
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
         B_p,K_p,P_sat_ratio] = dynamics_error_wrapper(X_f_guess,X_s_guess,...
                                                       phase_X_f_guess,phase_X_s_guess,...
                                                       ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                                drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody);

        X_err       = max( abs([X_f_err,X_s_err]),             [], 'all');
        phase_X_err = max( abs([phase_X_f_err,phase_X_s_err]), [], 'all');
        F_lim_err   = max( abs(force_lim_err),                 [], 'all');
        X_lim_err   = max( abs(amp_lim_err),                   [], 'all');

        % new guesses
        X_f_guess = mag_X_f;
        X_s_guess = mag_X_s;
        phase_X_f_guess = phase_X_f;
        phase_X_s_guess = phase_X_s;
        ctrl_mult_guess = ctrl_mult_guess ./ (force_lim_err+1);
        phase_ctrl_mult_guess = eps + zeros(size(ctrl_mult_guess));

        % check convergence
        X_converged = X_err < X_tol;
        phase_X_converged = phase_X_err < phase_X_tol || ...
                abs(phase_X_err - 2*pi) < phase_X_tol ; % the 2pi allows for the solution to oscillate between +pi and -pi
        F_lim_converged = F_lim_err < 0.01; %F_lim_tol;
        X_lim_converged = X_lim_err < 0.01; %X_lim_tol;
        converged = X_converged && phase_X_converged && ...
                    F_lim_converged && X_lim_converged;

        if all(isnan(X_f_guess))
            error('all nan')
        end

        % increment iterations and check max iters
        iters = iters + 1;
        if iters > max_drag_iters
%             warning(['Drag loop has not converged after ' num2str(max_drag_iters) ' iterations. ' ...
%                 'X_err = ' num2str(X_err) ', phase_X_err = ' num2str(phase_X_err)])
            break
        end
    end

    if drag_convergence_plot_on % fixme: have guesses as output of dynamics and save it to vals and plot outside of the sim
        plot_drag_convergence(X_f_guesses, X_s_guesses, phase_X_f_guesses, phase_X_s_guesses, ...
                              iters, multibody, X_tol, phase_X_tol)
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
          B_p,K_p,P_sat_ratio] = dynamics_error_wrapper(X_f_guess,X_s_guess,...
                                                        phase_X_f_guess,phase_X_s_guess,...
                                                        ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                                F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                                drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                                X_max,control_type,multibody)
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p,P_sat_ratio,...
     force_lim_err,amp_lim_err] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                    X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                    B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,w,...
                                    F_f_mag,F_f_phase,F_s_mag,F_s_phase,...
                                    control_type,multibody,F_max,X_max,...
                                    ctrl_mult_guess,phase_ctrl_mult_guess); % closed loop dynamics
    % error
    X_f_err = X_f_guess - mag_X_f;
    X_s_err = X_s_guess - mag_X_s;
    phase_X_f_err = phase_X_f_guess - phase_X_f;
    phase_X_s_err = phase_X_s_guess - phase_X_s;

end

function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p_sat,K_p_sat,...
         P_sat_ratio,...
         force_lim_err,...
         amp_lim_err] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                        X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                        B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,w,...
                                        F_f_mag,F_f_phase,F_s_mag,F_s_phase,...
                                        control_type,multibody,F_max,X_max,...
                                        ctrl_mult_guess,phase_ctrl_mult_guess)

    [B_drag_f, K_drag_f] = get_drag_dynamic_coeffs(X_f_guess, phase_X_f_guess, mag_v0_f, w, drag_const_f);
    [B_drag_s, K_drag_s] = get_drag_dynamic_coeffs(X_s_guess, phase_X_s_guess, mag_v0_s, w, drag_const_s);

    B_f = B_h_f + B_drag_f;
    B_s = B_h_s + B_drag_s;
    K_f = K_h_f + K_drag_f;
    K_s = K_h_s + K_drag_s;

    % unsaturated control gains
    [B_p,K_p] = controller(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w,control_type,multibody);

    % amplitude and power output for chosen controller
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p_sat,K_p_sat,...
     P_sat_ratio,...
     force_lim_err,...
     amp_lim_err] = get_response_saturated(B_c,B_f,B_s,K_f,K_s,...
                                            m_c,m_f,m_s,w,K_p,B_p,...
                                            F_f_mag,F_f_phase,...
                                            F_s_mag,F_s_phase,...
                                            F_max,X_max,multibody,...
                                            ctrl_mult_guess,phase_ctrl_mult_guess);
end

function [B_p,K_p] = controller(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w, control_type, multibody)

    % intrinsic admittance of the controlled DOF
    if multibody
        [real_G_u,imag_G_u] = multibody_impedance(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w);
        mag_G_u_squared = real_G_u.^2 + imag_G_u.^2;
    else
        % derived on page 58 of notebook
        resistance = B_f;
        reactance = m_f.*w - K_f./w;
        mag_G_u_squared = 1 ./ (resistance.^2 + reactance.^2);
        real_G_u = resistance .* mag_G_u_squared;
        imag_G_u = -reactance .* mag_G_u_squared;
    end

    % control - set powertrain coefficients
    if strcmpi(control_type,'reactive')
        K_p = -w .* imag_G_u ./ mag_G_u_squared;
        B_p = real_G_u ./ mag_G_u_squared;
    elseif strcmpi(control_type, 'damping')
        B_p = 1 ./ sqrt(mag_G_u_squared);
        K_p = 1e-8; % can't be quite zero because r_k = Inf
    end
    B_p = max(B_p,0); % don't allow negative damping
end

function [B_drag, K_drag] = get_drag_dynamic_coeffs(X_guess, phase_X_guess, mag_v0, w, drag_const)
    mag_v = w .* X_guess;
    phase_v = phase_X_guess + pi/2;
    mag_v0_v_ratio = mag_v0 ./ mag_v;
    mag_v0_v_ratio(mag_v==0) = 0; % prevent divide by zero
    phase_v_prime = atan2( cos(phase_X_guess) - mag_v0_v_ratio, -sin(phase_X_guess)); % derived on p67 of notebook #6 (10/4/24)
    
    alpha_v = sqrt(1 + mag_v0_v_ratio.^2 - 2 * mag_v0_v_ratio .* cos(phase_X_guess)); % derived on p48 of notebook #5 (6/7/24)
    phi_alpha = 0; %phase_v_prime - phase_v;
    mag_cf = drag_const * alpha_v.^2 .* mag_v; % eq 52 in Water paper

    B_drag = mag_cf .* cos(phi_alpha); % real part of c_f
    B_drag = max(B_drag,zeros(size(B_drag))); % prevent negative damping
    K_drag = - w .* mag_cf .* sin(phi_alpha); % -w times imag part of c_f
end

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p_sat,K_p_sat,...
     P_sat_ratio,...
     F_err,X_err] = get_response_saturated(B_c,B_f,B_s,K_f,K_s,...
                                                m_c,m_f,m_s,w,K_p,B_p,...
                                                F_f_mag,F_f_phase,...
                                                F_s_mag,F_s_phase,...
                                                F_max,X_max,multibody,...
                                                ctrl_mult_guess,...
                                                phase_ctrl_mult_guess)
    if multibody
        [mag_U_unsat,~,P_unsat,...
         ~,~,~,mag_X_f_unsat] = multibody_response(B_c, B_f, B_s, K_f, K_s,...
                                         m_c, m_f, m_s, w, K_p, B_p,...
                                         F_f_mag, F_f_phase,...
                                         F_s_mag, F_s_phase);
    else
        b = B_f + B_p;
        k = K_f + K_p;
        mag_X_f_unsat = second_order_transfer_fcn(w, m_f, b, k, F_f_mag);
        F_ptrain_over_x = sqrt( (B_p .* w).^2 + (K_p).^2 );
        mag_U_unsat = F_ptrain_over_x .* mag_X_f_unsat;
        P_unsat = 1/2 * B_p .* w.^2 .* mag_X_f_unsat.^2;
    end

    % get force-saturated response



    % fixme: If reactive control, mult should be complex, based on eq4 of IFAC paper.
%     F_err = zeros(size(f_sat));
%     max_err = 1;
%     mult = f_sat;
%     iters = 0;

%     err_thresh = 0.01;
%     while max_err > err_thresh
% 
%         iters = iters + 1;
%         if iters > 50
%             warning('force saturation loop failed to converge. reversing search direction on problem sea states.')
%             F_err(abs(F_err) > err_thresh) = -F_err(abs(F_err) > err_thresh);
%         end
%         if iters > 100
%             warning('force saturation loop failed to converge, and reversing search didnt help')
%             break
%         end
%         mult = mult ./ (F_err+1);%get_multiplier(f_sat,m_f,B_f,K_f,w, B_f./B_p, K_f./K_p); % fixme this is wrong for multibody
    
        mult = ctrl_mult_guess .* exp(1i*phase_ctrl_mult_guess);
        Z_p = B_p + K_p ./ (1i * w);
        Z_p_sat = Z_p .* mult;

        B_p_sat =       real(Z_p_sat);
        K_p_sat = -w .* imag(Z_p_sat);
    
        if multibody
            [mag_U,phase_U,...
             real_P,reactive_P,...
             mag_X_u,phase_X_u,...
             mag_X_f,phase_X_f,...
             mag_X_s,phase_X_s] = multibody_response(B_c, B_f, B_s, K_f, K_s, ...
                                                     m_c, m_f, m_s, w, ...
                                                     K_p_sat, B_p_sat, ...
                                                     F_f_mag, F_f_phase, ...
                                                     F_s_mag, F_s_phase);
        else
            b_sat = B_f + B_p_sat;
            k_sat = K_f + K_p_sat;
    
            [mag_X_f,phase_X_f] = second_order_transfer_fcn(w, m_f, b_sat, k_sat, F_f_mag, F_f_phase);
            mag_X_u = mag_X_f;   phase_X_u = phase_X_f;
            mag_X_s = 0*mag_X_f; phase_X_s = 0*phase_X_f;
            mag_U = ctrl_mult_guess .* F_ptrain_over_x .* mag_X_f;
            
            phase_Z_u = atan2(-K_p_sat./w, B_p_sat);    % phase of control impedance
            phase_V_u = pi/2 + phase_X_u;               % phase of control velocity
            phase_U = phase_V_u + phase_Z_u;            % phase of control force
    
            real_P = 1/2 * B_p_sat .* w.^2 .* mag_X_u.^2; % this is correct even if X and U are out of phase
            check_P = 1/2 * w .* mag_X_u .* mag_U .* cos(phase_U - phase_V_u); % so is this, they match
            reactive_P = 0; % fixme this is incorrect but doesn't affect anything rn
        end
    

        % notebook p106 2/2/25
%         f_sat = min(4/pi * F_max ./ mag_U_unsat,   1);
%         x_sat = min(       X_max ./ mag_X_f_unsat, 1);

        % what the force and amp would be if the force-sat and amp-sat
        % solutions applied
        force_saturated_U = min(4/pi * F_max, mag_U_unsat);
        amp_saturated_X   = min(       X_max, mag_X_f_unsat);

        % indices where sat and unsat solutions violate which constraint
%         idx_force_viol_sat   = mag_U         > force_saturated_U;
%         idx_force_viol_unsat = mag_U_unsat   > force_saturated_U;
%         idx_amp_viol_sat     = mag_X_f       > amp_saturated_X;
%         idx_amp_viol_unsat   = mag_X_f_unsat > amp_saturated_X;

%         both_ok_unsat       = ~idx_force_viol_unsat & ~idx_amp_viol_unsat;
%         only_force_ok_unsat = ~idx_force_viol_unsat &  idx_amp_viol_unsat;
%         only_amp_ok_unsat   =  idx_force_viol_unsat & ~idx_amp_viol_unsat;

        % indices where each solution applies. currently overriding amp
        % constraint. commented forumlas attempt to incorporate both but
        % should not actually be based purely on whether the 
        % saturated and/or unsaturated solution violates, it should also
        % depend on alpha and the limits (see notebook p142-144 9/15/25).
        if isinf(F_max)
            idx_force_sat_applies = false(size(w));
        else
            idx_force_sat_applies = true(size(w));
        end
        idx_amp_sat_applies   = false(size(w));
%         idx_force_sat_applies = both_ok_unsat | idx_force_viol_sat;   % | (idx_amp_viol_sat & idx_force_viol_sat);
%         idx_amp_sat_applies   = both_ok_unsat | only_force_ok_unsat; % | (idx_amp_viol_sat & idx_force_viol_sat);

        F_err_from_force_sat = mag_U   ./ force_saturated_U - 1;
        X_err_from_amp_sat   = mag_X_f ./ amp_saturated_X   - 1;

        F_err = Inf(size(F_err_from_force_sat));
        X_err = Inf(size(X_err_from_amp_sat));

        % at sea states where the force-saturated solution applies, 
        % set F_err as deviation from that solution. Otherwise, set F_err=0
        % to allow any force to occur.
        F_err(idx_force_sat_applies) = F_err_from_force_sat(idx_force_sat_applies);
        F_err(~idx_force_sat_applies) = 0;

        % at sea states where the amplitude-saturated solution applies, 
        % set X_err as deviation from that solution. Otherwise, set X_err = 0
        % to allow any amplitude to occur.
        X_err(idx_amp_sat_applies) = X_err_from_amp_sat(idx_amp_sat_applies);
        X_err(~idx_amp_sat_applies) = 0;

        % prevent 0/0=NaN when limit of zero is set
        if F_max==0
            F_err(mag_U==0) = 0;
            F_err(mag_U~=0) = mag_U(mag_U~=0) ./ mag_U_unsat(mag_U~=0);
        end

        % add penalty for B_p<0 (negative power)
        tol = 1;
        B_p_violation = max(-B_p_sat+tol,0);
        F_err = F_err + B_p_violation;

    P_sat_ratio = real_P ./ P_unsat;

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

function mult = get_multiplier(f_sat,m,b,k,w,r_b,r_k)
    % m, k, and r_k are scalars.
    % All other inputs are 2D arrays, the dimension of the sea state matrix.

    % speedup: only do math for saturated sea states, since unsat will = 1
    % likewise, don't do math for uncontrolled sea states (f_sat=0)
    idx_no_sat = f_sat == 1;
    idx_zero = f_sat == 0;
    idx_nan = idx_no_sat | idx_zero;
    f_sat(idx_nan) = NaN;
    b(idx_nan) = NaN;
    w(idx_nan) = NaN;
    r_b(idx_nan) = NaN;

    [a_quad, b_quad, c_quad]  = get_abc_symbolic(f_sat,m,b,k,w,r_b,r_k);

    % solve the quadratic formula
    determinant = sqrt(b_quad .^ 2 - 4 * a_quad .* c_quad);
    num = -b_quad + determinant;
    num(:,:,2) = -b_quad - determinant;
    den = 2 * a_quad;
    roots = num ./ den;

    % choose which of the two roots to use
    mult = pick_which_root(roots, idx_no_sat, idx_zero, a_quad, b_quad, c_quad);
    assert(all(~isnan(mult),'all'))
end