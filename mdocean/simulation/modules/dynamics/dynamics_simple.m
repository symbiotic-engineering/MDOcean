function [w,A,B,K,Fd,k] = dynamics_simple(Hs, T, D_f, T_f, D_s, T_s, h, rho_w, g, use_MEEM, harmonics)
    w = 2*pi./T;        % angular frequency
    k = w.^2 / g;       % wave number (dispersion relation for deep water)

    a2 = D_f / 2;        % radius
    a1 = D_s / 2;
    d2 = T_f;        % draft below waterline
    d1 = T_s;
    A_w = pi * (a2^2 - a1^2);     % waterplane area
    
    if use_MEEM
        [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs_MEEM(a2, k, d2, a1, d1, h, harmonics); 
    else
        [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs(a2, k, d2);  
    end

    A = rho_w * A_over_rho;                 % added mass
    B = rho_w * w .* B_over_rho_w;          % radiation damping
    gamma   = rho_w * g * gamma_over_rho_g; % froude krylov coefficient
    K       = rho_w * g * A_w;              % hydrostatic stiffness
    Fd      = gamma .* Hs;                  % excitation force of wave
end