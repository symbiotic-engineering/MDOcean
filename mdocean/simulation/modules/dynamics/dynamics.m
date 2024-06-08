function [F_heave_max, F_surge_max, F_ptrain_max, ...
    P_var, P_elec, P_matrix, h_s_extra, P_unsat] = dynamics(in,m_float,V_d,draft)

    % use probabilistic sea states for power
    [T,Hs] = meshgrid(in.T,in.Hs);
    [P_matrix,h_s_extra,P_unsat] = get_power_force(in,T,Hs,m_float,V_d,draft);
    
    % account for powertrain electrical losses
    P_matrix = P_matrix * in.eff_pto;
    
    % saturate maximum power
    P_matrix = min(P_matrix,in.power_max);
    
    % weight power across all sea states
    P_weighted = P_matrix .* in.JPD / 100;
    P_elec = sum(P_weighted(:)); 
    
    assert(isreal(P_elec))
    
    % use max sea states for structural forces and max amplitude
    [~,~,~,F_heave_max,F_surge_max,...
        F_ptrain_max] = get_power_force(in, ...
                                in.T_struct, in.Hs_struct, m_float, V_d, draft);
    
    % coefficient of variance (normalized standard deviation) of power
    P_var = std(P_matrix(:), in.JPD(:)) / P_elec;
    P_var = P_var * 100; % convert to percentage

end

function [P_matrix, h_s_extra, P_unsat, ...
          F_heave, F_surge, F_ptrain_max] = get_power_force(in,T,Hs, m_float,V_d, draft)

    % get dynamic coefficients
    [w,A,B_h,K_h,F_d,k_wvn,drag_const] = dynamics_simple(Hs, T, in.D_f, in.T_f, in.D_s, ...
                                            in.T_s, in.h, in.C_d_float, in.rho_w, in.g, ...
                                            in.use_MEEM, in.harmonics);
    m = m_float + A;

    % control - set powertrain coefficients
    if strcmp(in.control_type,'reactive')
        K_p = w.^2 .* m - K_h;
        B_p = B_h;
    elseif strcmp (in.control_type, 'constant impedance')
        K_p = in.w_n^2 * m - K_h; % fixme this is not actually constant since m varies
        B_p = in.B_p;
    elseif strcmp (in.control_type, 'damping')
        B_p = sqrt(B_h.^2 + ( (w.^2 .* m - K_h) ./ w).^2);
        K_p = 1e-8; % can't be quite zero because r_k = Inf
    end

    % get response: includes drag and force saturation
    [X, X_unsat, mult, F_ptrain] = get_response_drag(w,m,B_h,B_p,K_h,K_p,F_d,in.F_max,drag_const);

    % confirm unsaturated response doesn't exceed maximum capture width
    P_unsat = 1/2 * B_p .* w.^2 .* X_unsat.^2;
    
    % calculate power
    P_matrix = 1/2 * (mult .* B_p) .* w.^2 .* X.^2;
    
    X_max = max(X,[],'all');
    % extra height on spar after accommodating float displacement
    h_s_extra_up = (in.h_s - in.T_s - (in.h_f - in.T_f) - X_max) / in.h_s;
    h_s_extra_down = (in.T_s - in.T_f - X_max) / in.h_s;
    h_s_extra = [h_s_extra_up, h_s_extra_down];

    % calculate forces
    if nargout > 3
        % powertrain force
        F_ptrain_max = max(F_ptrain,[],'all');
        F_ptrain_max = min(F_ptrain_max, in.F_max);

        % heave force: includes powertrain force and D'Alembert force
        F_heave_fund = sqrt( (mult .* B_p .* w).^2 + (mult .* K_p - m_float .* w.^2).^2 ) .* X;
        F_heave = min(F_heave_fund, in.F_max + m_float * w.^2 .* X);

        % surge force
        F_surge = max(Hs,[],'all') * in.rho_w * in.g * V_d .* (1 - exp(-max(k_wvn,[],'all')*draft));
    end
end

function [X, X_unsat, mult, F_ptrain] = get_response_drag(w,m,B_h,B_p,K_h,K_p,F_d,F_max,drag_const)
    % initial guess: no drag
    X_guess = 0;
    X_err = 1;
    iters = 0;

    while X_err > 1e-2
        alpha_v = 1; % fixme
        phi_alpha = 0; % fixme
        mag_cf = drag_const * alpha_v.^2 .* w .* X_guess; % eq 52 in Water paper

        B_drag = mag_cf .* cos(phi_alpha); % real part of c_f
        K_drag = - w .* mag_cf * sin(phi_alpha); % -w times imag part of c_f

        B_not_p = B_h + B_drag;
        K_not_p = K_h + K_drag;
        [X, X_unsat, mult, F_ptrain] = get_response_saturated(w,m,B_not_p,B_p,K_not_p,K_p,F_d,F_max);
        X_err = max( abs(X_guess - X), [], 'all');
        X_guess = X;
        iters = iters + 1;
        disp(iters);
    end
end

function [X_sat,X_unsat,mult,F_ptrain] = get_response_saturated(w,m,B_not_p,B_p,K_not_p,K_p,F_d,F_max)
% "not p" refers to anything that is not the powertrain under the force limit (so hydro, drag, drivetrain, etc)
    b = B_not_p + B_p;
    k = K_not_p + K_p;
    X_unsat = second_order_transfer_fcn(w, m, b, k, F_d);
    F_ptrain_over_x = sqrt( (B_p .* w).^2 + (K_p).^2 );
    F_ptrain_unsat = F_ptrain_over_x .* X_unsat;
    
    % get saturated response
    r = min(F_max ./ F_ptrain_unsat, 1);%fcn2optimexpr(@min, in.F_max ./ F_ptrain_unsat, 1);
    alpha = 2/pi * ( 1./r .* asin(r) + sqrt(1 - r.^2) );
    f_sat = alpha .* r;
    mult = get_multiplier(f_sat,m,b,k,w, b./B_p, k./K_p);
    b_sat = B_not_p + mult .* B_p;
    k_sat = K_not_p + mult .* K_p;
    X_sat = second_order_transfer_fcn(w, m, b_sat, k_sat, F_d);

    F_ptrain = mult .* F_ptrain_over_x .* X_sat;
    
    F_err_1 = abs(F_ptrain ./ (F_max * alpha) - 1);
    F_err_2 = abs(F_ptrain ./ (f_sat .* F_ptrain_unsat) - 1);
    % 0.1 percent error
    if any(f_sat<1,'all')
        stuff = all(F_err_1(f_sat < 1) < 1e-3, 'all');
        if ~stuff
            stuff
        end
        assert(stuff);
    end
    assert(all(F_err_2 < 1e-3,'all'));
end

function X = second_order_transfer_fcn(w,m,b,k,F)
    imag_term = b .* w;
    real_term = k - m .* w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    %X_over_F_phase = atan2(imag_term,real_term);
    X = X_over_F_mag .* F;
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