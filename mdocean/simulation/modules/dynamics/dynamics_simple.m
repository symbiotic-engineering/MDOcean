function [w,A,B,K,Fd,k] = dynamics_simple(Hs, T, D_f, T_f, rho_w, g)
    w = 2*pi./T;        % angular frequency
    k = w.^2 / g;       % wave number (dispersion relation for deep water)

    r = D_f / 2;        % radius
    draft = T_f;        % draft below waterline
    A_w = pi * r^2;     % waterplane area
    
    [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs(r, k, draft);  
    
    A = rho_w * A_over_rho;                 % added mass
    B = rho_w * w .* B_over_rho_w;          % radiation damping
    gamma   = rho_w * g * gamma_over_rho_g; % froude krylov coefficient
    K       = rho_w * g * A_w;              % hydrostatic stiffness
    Fd      = gamma .* Hs;                  % excitation force of wave
end