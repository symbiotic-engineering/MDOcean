
function [F_heave, F_surge, F_ptrain, P_var, P_elec, P_matrix] = dynamics(x,p,m_float,V_d,draft)

% use probabilistic sea states for power
[Hs,T] = meshgrid(p.T,p.Hs);
P_matrix = get_power_force(x,p,T,Hs,m_float,V_d,draft);

% weight power across all sea states
P_weighted = P_matrix .* p.JPD / 100;
P_elec = sum(P_weighted(:)); 

% use max sea states for structures
[~,F_heave,F_surge,F_ptrain] = get_power_force(x, p, p.T_struct, p.Hs_struct, m_float, V_d, draft);

% coefficient of variance (normalized standard deviation) of power
P_var = std(P_matrix(:), p.JPD(:)) / P_elec;
P_var = P_var * 100; % convert to percentage

end

function [P_matrix, F_heave, F_surge, F_ptrain] = get_power_force(x,p,T,Hs, m_float,V_d, draft)
    % get unsaturated response
    [w,A,B,K,Fd,k_wvn] = dynamics_simple(Hs, T, x.D_sft, p.rho_w, p.g);       
    m = m_float + A;
    b = B + x.D_int;
    k = x.w_n^2 * m;
    K_int = k - K;
    X_unsat = get_response(w,m,b,k,Fd);
    F_ptrain_unsat = sqrt( (x.D_int * w).^2 + (K_int).^2 ).* X_unsat;
    
    % get saturated response
    mult = min(x.F_max ./ F_ptrain_unsat, 1);%fcn2optimexpr(@min, x.F_max ./ F_ptrain_unsat, 1);
    % fixme: should multiply mult (saturation multiplier) by a fourier multiplier to get total mult
    b_sat = B + mult * x.D_int;
    k_sat = K + mult * K_int;
    X_sat = get_response(w,m,b_sat,k_sat,Fd);
    
    P_matrix = 1/2 * (mult * x.D_int) .* w.^2 .* X_sat.^2;
    
    if nargout > 1
        F_ptrain = mult .* sqrt( (x.D_int*w).^2 + K_int^2 )* X_sat; % todo: check that this doesn't exceed F_max
        F_heave = ((p.rho_w*p.g*(p.Hs/2)*pi)/4)*(((x.D_or^2*exp(-k_wvn*(draft(1)+draft(2)+p.t_r))))-((x.D_or^2-x.D_i^2)*exp(-k_wvn*(draft(1)+draft(2))))+((x.D_sft^2-x.D_i^2)*exp(-k_wvn*draft(1)))); added mass and excitation
        F_surge = Hs * p.rho_w * p.g * V_d * (1 - exp(-k_wvn*draft));
    end
end

function [w,A,B,K,Fd,k] = dynamics_simple(Hs, T, D_sft, rho_w, g)
    w = 2*pi./T;         % frequency
    k = w.^2 / g;        % wave number
    V_g = g ./(2*w);     % group velocity

    r = D_sft / 2;      % radius
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
