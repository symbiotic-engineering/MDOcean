function [F_heave_max, F_surge_max, F_ptrain_max, ...
    P_var, P_elec, P_matrix, h_s_extra, P_unsat, X_below_wave, X_below_linear] = ...
    dynamics(in,m_float,V_d,draft)

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
        F_ptrain_max,X_sat] = get_power_force(in, ...
                                in.T_struct, in.Hs_struct, m_float, V_d, draft);
    
    % coefficient of variance (normalized standard deviation) of power
    P_var = std(P_matrix(:), in.JPD(:)) / P_elec;
    P_var = P_var * 100; % convert to percentage
    
    % prevent rising out of the water
    X_max_linear = 1/10 * in.D_f;
    X_below_wave = Hs - X_sat;
    X_below_linear = X_max_linear - X_sat;
end

function [P_matrix, h_s_extra, P_unsat, F_heave, F_surge, F_ptrain_max, X_sat] = ...
get_power_force(in,T,Hs, m_float,V_d, draft)
    % get unsaturated response
    [w,A,B_h,K_h,Fd,k_wvn] = dynamics_simple(Hs, T, in.D_f, in.T_f, in.rho_w, in.g);
    m = m_float + A;
    if in.freq_based_optimal_ctrl
        in.B_p = B_h;
        k = (2*pi./T).^2 * m;
    else
        k = in.w_n^2 * m;
    end
    b = B_h + in.B_p;
    K_p = k - K_h;
    X_unsat = get_response(w,m,b,k,Fd);

    % confirm unsaturated response doesn't exceed maximum capture width
    P_unsat = 1/2 * in.B_p .* w.^2 .* X_unsat.^2;

    F_ptrain_over_x = sqrt( (in.B_p .* w).^2 + (K_p).^2 );
    F_ptrain_unsat = F_ptrain_over_x .* X_unsat;
    
    % get saturated response
    r = min(in.F_max ./ F_ptrain_unsat, 1);%fcn2optimexpr(@min, in.F_max ./ F_ptrain_unsat, 1);
    alpha = 2/pi * ( 1./r .* asin(r) + sqrt(1 - r.^2) );
    f_sat = alpha .* r;
    mult = get_multiplier(f_sat,m,b,k,w, b./in.B_p, k./K_p);
    b_sat = B_h + mult .* in.B_p;
    k_sat = K_h + mult .* K_p;
    X_sat = get_response(w,m,b_sat,k_sat,Fd);
    
    % calculate power
    P_matrix = 1/2 * (mult .* in.B_p) .* w.^2 .* X_sat.^2;
    X_max = max(X_sat,[],'all');
    h_s_extra = (in.h_s - in.T_s - (in.h_f - in.T_f) - X_max) / in.h_s; % extra height on spar after accommodating float displacement

    % calculate forces
    if nargout > 2
        F_ptrain = mult .* F_ptrain_over_x .* X_sat;
        F_ptrain_max = max(F_ptrain,[],'all');
        F_err_1 = abs(F_ptrain ./ (in.F_max * alpha) - 1);
        F_err_2 = abs(F_ptrain ./ (f_sat .* F_ptrain_unsat) - 1);
        % 0.1 percent error
        if any(f_sat<1,'all')
%            assert(all(F_err_1(f_sat < 1) < 1e-3),'all');
        end
%        assert(all(F_err_2 < 1e-3,'all'));

        F_heave_fund = sqrt( (mult .* in.B_p .* w).^2 + (mult .* K_p - m_float * w.^2).^2 ) .* X_sat; % includes powertrain force and D'Alembert force
        F_heave = min(F_heave_fund, in.F_max + m_float * w.^2 .* X_sat);
        %assert(F_heave <= in.F_max);

        F_surge = max(Hs,[],'all') * in.rho_w * in.g * V_d .* (1 - exp(-max(k_wvn,[],'all')*draft));
    end
end

function X = get_response(w,m,b,k,Fd)
    imag_term = b.*w;
    real_term = k - m*w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    %X_over_F_phase = atan2(imag_term,real_term);
    X = X_over_F_mag .* Fd;
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