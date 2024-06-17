function [m_f,B_h_f,K_h_f,F_f_mag,F_f_phase,...
         m_s,B_h_s,K_h_s,F_s_mag,F_s_phase,...
         m_c,B_c,drag_const_f,drag_const_s,...
         mag_v0_f,mag_v0_s,w,k] = get_dynamic_coeffs(Hs, T, ...
                                                D_f, T_f, D_s, D_d, T_s, h, ...
                                                m_float, m_spar, spar_excitation_coeffs,...
                                                C_d_float, C_d_spar, ...
                                                rho_w, g, mu, ...
                                                use_MEEM, harmonics, hydro)
    w = 2*pi./T;        % angular frequency
    k = w.^2 / g;       % wave number (dispersion relation for deep water)

    a2 = D_f / 2;        % radius
    a1 = D_s / 2;
    d2 = T_f;        % draft below waterline
    d1 = T_s;
    A_w_f = pi * (a2^2 - a1^2);     % waterplane area
    
    if use_MEEM
        [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs_MEEM(a2, k, d2, a1, d1, h, harmonics); 
    else
        [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs(a2, k, d2, hydro);  
    end

    A_f = rho_w * A_over_rho;               % added mass
    m_f = A_f + m_float;
    B_h_f   = rho_w * w .* B_over_rho_w;    % radiation damping
    gamma   = rho_w * g * gamma_over_rho_g; % froude krylov coefficient
    K_h_f   = rho_w * g * A_w_f;              % hydrostatic stiffness
    H       = Hs / sqrt(2);                 % equivalent regular wave height
    F_f_mag = gamma .* H / 2;               % excitation force of wave

    F_f_phase = 0; F_s_phase = 0; % fixme excitation phase
    m_c = 0; B_c = 0;             % fixme coupling

    [m_s,K_h_s,F_s_mag,B_h_s] = spar_dynamics(D_s/D_d, D_d, rho_w, mu, C_d_spar, T_s, ...
                                                spar_excitation_coeffs, m_spar, Hs, w, g);

    A_drag_s = pi/4 * D_d^2;

    drag_const_f = 4/(3*pi) * rho_w * A_w_f * C_d_float;
    drag_const_s = 4/(3*pi) * rho_w * A_drag_s * C_d_spar; % fixme use KC as in spar_dynamics
    mag_v0_f = H/2 * g .* k ./ w .* exp(-k * T_f/2); % finite depth velocity of incident wave
    mag_v0_s = H/2 * g .* k ./ w .* exp(-k * T_s/2);
end