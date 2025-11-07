% called by r_k_zeta - has been replaced by get_dynamic_coeffs  / get_dynamic_coeffs_MEEM
function [w, A, B, K, Fd, k, drag_const, mag_v0] = dynamics_simple(Hs, T, D_f, T_f, D_s, T_s, h, C_d, rho_w, g, use_MEEM, harmonics, hydro)
    w = 2 * pi ./ T;        % angular frequency
    k = w.^2 / g;       % wave number (dispersion relation for deep water)

    a2 = D_f / 2;        % radius
    a1 = D_s / 2;
    d2 = T_f;        % draft below waterline
    d1 = T_s;
    A_w = pi * (a2^2 - a1^2);     % waterplane area

    if use_MEEM
        [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs_MEEM(a2, k, d2, a1, d1, h, harmonics);
    else
        [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs(a2, k, d2, hydro);
    end

    A = rho_w * A_over_rho;                 % added mass
    B = rho_w * w .* B_over_rho_w;          % radiation damping
    gamma   = rho_w * g * gamma_over_rho_g; % froude krylov coefficient
    K       = rho_w * g * A_w;              % hydrostatic stiffness
    H       = Hs / sqrt(2);                 % equivalent regular wave height
    Fd      = gamma .* H / 2;               % excitation force of wave

    drag_const = 4 / (3 * pi) * rho_w * A_w * C_d;
    mag_v0 = H / 2 * g .* k ./ w .* exp(-k * T_f / 2); % finite depth velocity of incident wave
end
