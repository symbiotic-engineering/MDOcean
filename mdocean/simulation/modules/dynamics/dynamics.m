function [F_heave_max, F_surge_max, F_ptrain_max, ...
    P_var, P_avg_elec, P_matrix_elec, X_constraints, B_p, X_u, P_matrix_mech] = dynamics(in,m_float,m_spar,V_d,draft)

    % use probabilistic sea states for power and PTO force and max amplitude
    [T,Hs] = meshgrid(in.T,in.Hs);
    [P_matrix_mech,X_constraints,B_p,X_u,~,~,F_ptrain_max] = get_power_force(in,T,Hs,m_float,m_spar,V_d,draft);
    
    % account for powertrain electrical losses
    P_matrix_elec = P_matrix_mech * in.eff_pto;
    
    % saturate maximum power
    P_matrix_elec = min(P_matrix_elec,in.power_max);
    
    % weight power across all sea states
    P_weighted = P_matrix_elec .* in.JPD / 100;
    P_avg_elec = sum(P_weighted(:)); 
    
    assert(isreal(P_avg_elec))
    
    % use max sea states for structural forces
    [~,~,~,~,F_heave_max,F_surge_max,~] = get_power_force(in, ...
                                in.T_struct, in.Hs_struct, m_float, m_spar, V_d, draft);
    
    % coefficient of variance (normalized standard deviation) of power
    P_var = std(P_matrix_elec(:), in.JPD(:)) / P_avg_elec;
    P_var = P_var * 100; % convert to percentage
    

end

function [P_matrix, X_constraints, B_p, mag_X_u,...
          F_heave_f, F_surge, F_ptrain_max] = get_power_force(in,T,Hs, m_float, m_spar, V_d, draft)

    % get dynamic coefficients for float and spar
    [m_f,B_h_f,K_h_f,F_f_mag,F_f_phase,...
     m_s,B_h_s,K_h_s,F_s_mag,F_s_phase,...
     m_c,B_c,drag_const_f,drag_const_s,...
     mag_v0_f,mag_v0_s,w,k_wvn] = get_dynamic_coeffs(Hs, T, ...
                                            in.D_f, in.T_f_2, in.D_s, in.D_d, in.T_s, in.h, ...
                                            m_float, m_spar, in.spar_excitation_coeffs,...
                                            in.C_d_float, in.C_d_spar, ...
                                            in.rho_w, in.g, ...
                                            in.use_MEEM, in.harmonics, in.hydro);

    X_max = 1e6;%min(Hs / (2*sqrt(2)), in.T_f);

    % get response: includes drag and force saturation
    [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p,K_p] = get_response_drag(w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                            F_f_mag,F_f_phase,F_s_mag,F_s_phase,in.F_max,...
                                            drag_const_f,drag_const_s,mag_v0_f,mag_v0_s, ...
                                            X_max,in.control_type,in.use_multibody,...
                                            in.X_tol,in.phase_X_tol,in.max_drag_iters);

% FIXME: check stability of closed loop multibody system
    
    % calculate power
    P_matrix = real_P;
    
    X_max = max(mag_X_u,[],'all');
    % extra height on spar after accommodating float displacement
    h_s_extra_up = (in.h_s - in.T_s - (in.h_f - in.T_f_2) - X_max) / in.h_s;
    h_s_extra_down = (in.T_s - in.T_f_2 - X_max) / in.h_s;

    % prevent rising out of the water
    X_max_linear = 1/10 * in.D_f;
    wave_amp = Hs(end,:)/(2*sqrt(2));
    X_below_wave = wave_amp ./ mag_X_u(end,:) - 1;
    X_below_linear = X_max_linear / X_max - 1;

    X_constraints = [h_s_extra_up, h_s_extra_down, X_below_linear, X_below_wave];

    % calculate forces
    if nargout > 3
        % powertrain force
        F_ptrain_max = max(mag_U,[],'all');
        F_ptrain_max = min(F_ptrain_max, in.F_max);

        % heave force: includes powertrain force and D'Alembert force
        %F_heave_fund = sqrt( (mult .* B_p .* w).^2 + (mult .* K_p - m_float .* w.^2).^2 ) .* X;
        
        F_heave_f = combine_ptrain_dalembert_forces(m_float, w, mag_X_f, phase_X_f, mag_U, phase_U, in.F_max);
        F_heave_s = combine_ptrain_dalembert_forces(m_spar,  w, mag_X_s, phase_X_s, mag_U, phase_U, in.F_max);

        % surge force
        F_surge = max(Hs,[],'all') * in.rho_w * in.g * V_d .* (1 - exp(-max(k_wvn,[],'all')*draft));
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
         B_p,K_p] = get_response_drag(w,m_f,m_s,m_c,B_h_f,B_h_s,B_c,K_h_f,K_h_s,...
                                        F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
                                        drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
                                        X_max,control_type,multibody,...
                                        X_tol,phase_X_tol,max_drag_iters)
    % initial guess: 2m float amplitude, 0.5m spar amplitude
    X_f_guess = 2;
    phase_X_f_guess = 0;
    X_s_guess = .5;
    phase_X_s_guess = 0;

    [X_err,phase_X_err] = deal(1);
    iters = 0;

    % loop until converged
    while X_err > X_tol || phase_X_err > phase_X_tol
        [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p,K_p] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
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

        % check iterations
        iters = iters + 1;
        if iters > max_drag_iters
            warning(['Drag loop has not converged after ' num2str(max_drag_iters) ' iterations. ' ...
                'X_err = ' num2str(X_err) ', phase_X_err = ' num2str(phase_X_err)])
            break
        end
    end
end

function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p_sat,K_p_sat] = dynamics_from_guess(X_f_guess, phase_X_f_guess, mag_v0_f, drag_const_f, ...
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
     B_p_sat,K_p_sat] = get_response_saturated(B_c,B_f,B_s,K_f,K_s,...
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
    phase_v_prime = atan2( cos(phase_X_guess) - mag_v0_v_ratio, -sin(phase_X_guess)); % derived on p67 of my notebook
    
    alpha_v = sqrt(1 + mag_v0_v_ratio.^2 - 2 * mag_v0_v_ratio .* cos(phase_X_guess)); % derived on p48 of my notebook
    phi_alpha = phase_v_prime - phase_v;
    mag_cf = drag_const * alpha_v.^2 .* mag_v; % eq 52 in Water paper

    B_drag = mag_cf .* cos(phi_alpha); % real part of c_f
    K_drag = - w .* mag_cf .* sin(phi_alpha); % -w times imag part of c_f
end

function [mag_U,phase_U,...
     real_P,reactive_P,...
     mag_X_u,phase_X_u,...
     mag_X_f,phase_X_f,...
     mag_X_s,phase_X_s,...
     B_p_sat,K_p_sat] = get_response_saturated(B_c,B_f,B_s,K_f,K_s,...
                                                m_c,m_f,m_s,w,K_p,B_p,...
                                                F_f_mag,F_f_phase,...
                                                F_s_mag,F_s_phase,...
                                                F_max,X_max,multibody)
    if multibody
        mag_U_unsat = multibody_response(B_c, B_f, B_s, K_f, K_s,...
                                         m_c, m_f, m_s, w, K_p, B_p,...
                                         F_f_mag, F_f_phase,...
                                         F_s_mag, F_s_phase);
    else
        b = B_f + B_p;
        k = K_f + K_p;
        X_unsat = second_order_transfer_fcn(w, m_f, b, k, F_f_mag);
        F_ptrain_over_x = sqrt( (B_p .* w).^2 + (K_p).^2 );
        mag_U_unsat = F_ptrain_over_x .* X_unsat;
    end

    % get force-saturated response
    r = min(F_max ./ mag_U_unsat, 1);%fcn2optimexpr(@min, in.F_max ./ F_ptrain_unsat, 1);
    alpha = 2/pi * ( 1./r .* asin(r) + sqrt(1 - r.^2) );
    f_sat = alpha .* r;
    mult = get_multiplier(f_sat,m_f,B_f,K_f,w, B_f./B_p, K_f./K_p); % fixme this is wrong for multibody

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

    F_err_1 = abs(mag_U ./ (F_max * alpha) - 1);
    F_err_2 = abs(mag_U ./ (f_sat .* mag_U_unsat) - 1);
    % 0.1 percent error
%     if any(f_sat<1,'all')
%         stuff = all(F_err_1(f_sat < 1) < 1e-3, 'all');
%         if ~stuff
%             stuff
%         end
%         assert(stuff);
%     end
%     assert(all(F_err_2 < 1e-3,'all'));

    % saturate position - hacky for now
    r_x = min(X_max ./ mag_X_f, 1);
    mag_X_u = mag_X_u .* r_x;
    mag_X_f = mag_X_f .* r_x;
    mag_X_s = mag_X_s .* r_x;
    mag_U   = mag_U   .* r_x;
    real_P  = real_P  .* r_x.^2;
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
    idx_no_sat = f_sat == 1;
    f_sat(idx_no_sat) = NaN;
    b(idx_no_sat) = NaN;
    w(idx_no_sat) = NaN;
    r_b(idx_no_sat) = NaN;
    
    [a_quad, b_quad, c_quad]  = get_abc_symbolic(f_sat,m,b,k,w,r_b,r_k);

    % solve the quadratic formula
    determinant = sqrt(b_quad .^ 2 - 4 * a_quad .* c_quad);
    num = -b_quad + determinant;
    num(:,:,2) = -b_quad - determinant;
    den = 2 * a_quad;
    roots = num ./ den;

    % choose which of the two roots to use
    mult = pick_which_root(roots, idx_no_sat, a_quad, b_quad, c_quad);
    assert(all(~isnan(mult),'all'))
end