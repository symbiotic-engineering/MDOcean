
function [F_heave_max, F_surge_max, F_ptrain_max, P_var, P_elec, P_matrix, h_s_extra] = dynamics(in,m_float,V_d,draft)

% use probabilistic sea states for power
[Hs,T] = meshgrid(in.T,in.Hs);
P_matrix = get_power_force(in,T,Hs,m_float,V_d,draft);

% weight power across all sea states
P_weighted = P_matrix .* in.JPD / 100;
P_elec = sum(P_weighted(:)); 

% use max sea states for structural forces and max amplitude
[~,F_heave_max,F_surge_max,...
    F_ptrain_max,h_s_extra] = get_power_force(in, ...
                            in.T_struct, in.Hs_struct, m_float, V_d, draft);

% coefficient of variance (normalized standard deviation) of power
P_var = std(P_matrix(:), in.JPD(:)) / P_elec;
P_var = P_var * 100; % convert to percentage

end

function [P_matrix, F_heave, F_surge, F_ptrain, F_ptrain_max, h_s_extra] = get_power_force(in,T,Hs, m_float,V_d, draft)
    % get unsaturated response
    [w,A,B_h,K_h,Fd,k_wvn] = dynamics_simple(Hs, T, in.D_f, in.rho_w, in.g);       
    m = m_float + A;
    b = B_h + in.B_p;
    k = in.w_n^2 * m;
    K_p = k - K_h;
    X_unsat = get_response(w,m,b,k,Fd);
    F_ptrain_over_x = sqrt( (in.B_p * w).^2 + (K_p).^2 );
    F_ptrain_unsat = F_ptrain_over_x .* X_unsat;
    
    % get saturated response
    r = in.F_max ./ F_ptrain_unsat;
    alpha = 2/pi * ( 1./r .* asin(r) + sqrt(1 - r.^2) );
    f_sat = min(alpha .* r, 1);%fcn2optimexpr(@min, in.F_max ./ F_ptrain_unsat, 1);
    mult = get_multiplier(f_sat,m,b,k,w);
    b_sat = B_h + mult * in.B_p;
    k_sat = K_h + mult * K_p;
    X_sat = get_response(w,m,b_sat,k_sat,Fd);
    
    P_matrix = 1/2 * (mult * in.B_p) .* w.^2 .* X_sat.^2;
    
    if nargout > 1
        F_ptrain = mult .* F_ptrain_over_x .* X_sat;
        F_ptrain_max = max(F_ptrain,[],'all');
        assert(F_ptrain == in.F_max * alpha);
        assert(F_ptrain == f_sat * F_ptrain_unsat);

        F_heave_fund = sqrt( (mult * in.B_p * w).^2 + (mult * K_p - m_float * w.^2).^2 ) .* X_sat; % includes powertrain force and D'Alembert force
        F_heave = min(F_heave_fund, in.F_max + m_float * w.^2 .* X_sat);
        assert(F_heave <= in.F_max);

        F_surge = Hs * in.rho_w * in.g * V_d .* (1 - exp(-k_wvn*draft));
        X_max = max(X_sat,[],'all');
        h_s_extra = (in.h_s - in.T_s - (in.h_f - in.T_f) - X_max) / in.h_s; % extra height on spar after accommodating float displacement
    end
end

function [w,A,B,K,Fd,k] = dynamics_simple(Hs, T, D_f, rho_w, g)
    w = 2*pi./T;         % frequency
    k = w.^2 / g;        % wave number
    V_g = g ./(2*w);     % group velocity

    r = D_f / 2;      % radius
    A_w = pi * r^2;     % waterplane area
    
    A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63; % added mass
    gamma   = rho_w * g * A_w; % Froude Krylov / diffraction
    B       = k ./ (4 * rho_w * g * V_g) * gamma.^2; % radiation damping
    K       = rho_w * g * A_w;  % hydrostatic stiffness
    Fd      = gamma * Hs;       % excitation force of wave
end

function X = get_response(w,m,b,k,Fd)
    real_term = b.*w;
    imag_term = k - m*w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    %X_over_F_phase = atan2(imag_term,real_term);
    X = X_over_F_mag .* Fd;
end

function mult = get_multiplier(f_sat,m,b,k,w)
    % all inputs are 2D arrays, the dimension of the sea state matrix
    
    % quadratic formula coeffs: 0 = a_quad * mult^2 + b_quad * mult + c_quad
    % this algebra is written out on p196 of my notebook
    m2_w4 = m.^2 .* w.^4;
    two_k_m_w2 = 2 * k .* m .* w.^2;
    f_sat2 = f_sat.^2;
    a_quad = (f_sat2 - 1) .* ((b.*w).^2 + k.^2) + two_k_m_w2 - m2_w4;
    b_quad = f_sat2 .* -two_k_m_w2;
    c_quad = f_sat2 .* m2_w4;

    % solve the quadratic formula
    determinant = sqrt(b_quad .^ 2 - 4 * a_quad .* c_quad);
    num = -b_quad + determinant;
    num(:,:,2) = -b_quad - determinant;
    den = 2 * a_quad;
    roots = num ./ den;

    % choose which of the two roots to use
    which_soln = roots == real(roots) & roots > 0 & roots <= 1.0001; % real solns on (0, 1]
    %both_ok = sum(which_soln,3) == 2;
    %which_soln(both_ok,2) = false; % if both are ok, choose the first one arbitrarily for now
    assert(all( sum(which_soln,3) == 1,'all') ); % confirm that 1 soln per sea state meets criteria
    
    mult = roots(which_soln);
    mult = reshape(mult,size(f_sat));

    % sanity check for the no saturation case
    err_no_sat = abs(mult(f_sat == 1) - 1);
    assert( all(err_no_sat < 1e-4,'all') ); 
end
