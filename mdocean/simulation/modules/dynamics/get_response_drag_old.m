function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio] = get_response_drag_old(w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                        F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                        drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                        X_max,control_type,multibody,...
                                        X_tol,phase_X_tol,max_drag_iters,drag_convergence_plot_on)
    % initial guess: 2m float amplitude, 0.5m spar amplitude
    X_f_guess = 2 * ones(size(w));
    phase_X_f_guess = zeros(size(w));
    X_s_guess = .5 * ones(size(w));
    phase_X_s_guess = zeros(size(w));

    converged = false;
    iters = 0;

    % loop until converged
    while ~converged
        if drag_convergence_plot_on
            % record guesses
            % 10,1 was found manually to be the problem sea state
            X_f_guesses(iters+1) = X_f_guess(10,1);
            X_s_guesses(iters+1) = X_s_guess(10,1);
            phase_X_f_guesses(iters+1) = phase_X_f_guess(10,1);
            phase_X_s_guesses(iters+1) = phase_X_s_guess(10,1);
        end

        [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p,P_sat_ratio] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                        X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                        B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,w,...
                                        F_f_mag,F_f_phase,F_s_mag,F_s_phase,...
                                        control_type,multibody,F_max,X_max); % closed loop dynamics
        % error
        X_f_err = X_f_guess - mag_X_f;
        X_s_err = X_s_guess - mag_X_s;
        phase_X_f_err = phase_X_f_guess - phase_X_f;
        phase_X_s_err = phase_X_s_guess - phase_X_s;
        X_err       = max( abs([X_f_err,X_s_err]),             [], 'all');
        phase_X_err = max( abs([phase_X_f_err,phase_X_s_err]), [], 'all');

        % new guesses
        X_f_guess = mag_X_f;
        X_s_guess = mag_X_s;
        phase_X_f_guess = phase_X_f;
        phase_X_s_guess = phase_X_s;

        % check convergence
        X_converged = X_err < X_tol;
        phase_X_converged = phase_X_err < phase_X_tol || phase_X_err < 2*pi+phase_X_tol; % the 2pi allows for the solution to oscillate between +pi and -pi
        converged = X_converged && phase_X_converged;

        % increment iterations and check max iters
        iters = iters + 1;
        if iters > max_drag_iters
            warning(['Drag loop has not converged after ' num2str(max_drag_iters) ' iterations. ' ...
                'X_err = ' num2str(X_err) ', phase_X_err = ' num2str(phase_X_err)])
            break
        end
    end

    if drag_convergence_plot_on % fixme: have guesses as output of dynamics and save it to vals and plot outside of the sim
        plot_drag_convergence(X_f_guesses, X_s_guesses, phase_X_f_guesses, phase_X_s_guesses, ...
                              iters, multibody, X_tol, phase_X_tol)
    end

end

function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p_sat,K_p_sat,...
         P_sat_ratio] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
                                        X_s_guess, phase_X_s_guess, mag_v0_s, drag_const_s, ...
                                        B_c,B_h_f,B_h_s,K_h_f,K_h_s,m_c,m_f,m_s,w,...
                                        F_f_mag,F_f_phase,F_s_mag,F_s_phase,...
                                        control_type,multibody,F_max,X_max)
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
     P_sat_ratio] = get_response_saturated(B_c,B_f,B_s,K_f,K_s,...
                                                m_c,m_f,m_s,w,K_p,B_p,...
                                                F_f_mag,F_f_phase,...
                                                F_s_mag,F_s_phase,...
                                                F_max,X_max,multibody);
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
     P_sat_ratio] = get_response_saturated(B_c,B_f,B_s,K_f,K_s,...
                                                m_c,m_f,m_s,w,K_p,B_p,...
                                                F_f_mag,F_f_phase,...
                                                F_s_mag,F_s_phase,...
                                                F_max,X_max,multibody)
    if multibody
        [mag_U_unsat,~,P_unsat] = multibody_response(B_c, B_f, B_s, K_f, K_s,...
                                         m_c, m_f, m_s, w, K_p, B_p,...
                                         F_f_mag, F_f_phase,...
                                         F_s_mag, F_s_phase);
    else
        b = B_f + B_p;
        k = K_f + K_p;
        X_unsat = second_order_transfer_fcn(w, m_f, b, k, F_f_mag);
        F_ptrain_over_x = sqrt( (B_p .* w).^2 + (K_p).^2 );
        mag_U_unsat = F_ptrain_over_x .* X_unsat;
        P_unsat = 1/2 * B_p .* w.^2 .* X_unsat.^2;
    end

    % get force-saturated response

    % notebook p106 2/2/25
    f_sat = min(4/pi * F_max ./ mag_U_unsat, 1);

    % fixme: If reactive control, mult should be complex, based on eq4 of IFAC paper.
    F_err = zeros(size(f_sat));
    max_err = 1;
    mult = f_sat;
    iters = 0;

    while max_err > 0.01

        iters = iters + 1;
        if iters > 100
            warning('force saturation loop failed to converge. reversing search direction on problem sea states.')
            F_err(abs(F_err) > max_err) = -F_err(abs(F_err) > max_err);
        end
        if iters > 1000
            warning('force saturation loop failed to converge, and reversing search didnt help')
            break
        end
        mult = mult ./ (F_err+1);%get_multiplier(f_sat,m_f,B_f,K_f,w, B_f./B_p, K_f./K_p); % fixme this is wrong for multibody
    
        B_p_sat = mult.*B_p;
        K_p_sat = mult.*K_p;
    
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
            mag_U = mult .* F_ptrain_over_x .* mag_X_f;
            
            phase_Z_u = atan2(-K_p_sat./w, B_p_sat);    % phase of control impedance
            phase_V_u = pi/2 + phase_X_u;               % phase of control velocity
            phase_U = phase_V_u + phase_Z_u;            % phase of control force
    
            real_P = 1/2 * B_p_sat .* w.^2 .* mag_X_u.^2; % this is correct even if X and U are out of phase
            check_P = 1/2 * w .* mag_X_u .* mag_U .* cos(phase_U - phase_V_u); % so is this, they match
            reactive_P = 0; % fixme this is incorrect but doesn't affect anything rn
        end
    
        F_err = mag_U ./ (f_sat .* mag_U_unsat) - 1;
        max_err = max(abs(F_err),[],'all');
    end

    % saturate position - hacky for now
    r_x = min(X_max ./ mag_X_f, 1);
    mag_X_u = mag_X_u .* r_x;
    mag_X_f = mag_X_f .* r_x;
    mag_X_s = mag_X_s .* r_x;
    mag_U   = mag_U   .* r_x;
    real_P  = real_P  .* r_x.^2;

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