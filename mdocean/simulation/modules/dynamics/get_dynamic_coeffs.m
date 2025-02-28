function [m_f,B_h_f,K_h_f,F_f_mag,F_f_phase,...
         m_s,B_h_s,K_h_s,F_s_mag,F_s_phase,...
         m_c,B_h_c,drag_const_f,drag_const_s,...
         mag_v0_f,mag_v0_s,w,k] = get_dynamic_coeffs(Hs, T, ...
                                                D_f, T_f, D_s, D_d, T_s, h, ...
                                                m_float, m_spar, spar_excitation_coeffs,...
                                                C_d_float, C_d_spar, ...
                                                rho_w, g, ...
                                                use_MEEM, harmonics, hydro)
    w = 2*pi./T;            % angular frequency
    k = dispersion(w,h,g);  % wave number (dispersion relation for finite depth water)

    a2 = D_f / 2;        % radius
    a1 = D_s / 2;
    a3 = D_d / 2;

    d2 = T_f;        % draft below waterline
    d1 = T_s;
    A_w_f = pi * (a2^2 - a1^2);     % waterplane area
    A_w_s = pi/4 * D_s^2; % fixme: for drag, should the area be the damping plate area, not the waterplane area?
    
    if use_MEEM
        [A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
        B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
        gamma_f_over_rho_g, gamma_s_over_rho_g, ...
        gamma_phase_f, gamma_phase_s] = get_hydro_coeffs_MEEM(a2, k, d2, a1, d1, a3, h, g, w, harmonics, spar_excitation_coeffs); 
    else
        [A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
        B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
        gamma_f_over_rho_g, gamma_s_over_rho_g, ...
        gamma_phase_f, gamma_phase_s] = get_hydro_coeffs(a2, k, d2, hydro);
    end

    A_f = rho_w * A_f_over_rho;               % added mass
    A_s = rho_w * A_s_over_rho;
    A_c = rho_w * A_c_over_rho;

    m_f = A_f + m_float;                      % total mass
    m_s = A_s + m_spar;
    m_c = A_c;

    B_h_f   = rho_w * w .* B_f_over_rho_w;    % radiation damping
    B_h_s   = rho_w * w .* B_s_over_rho_w;
    B_h_c   = rho_w * w .* B_c_over_rho_w;

    gamma_f = rho_w * g * gamma_f_over_rho_g; % excitation coefficient
    gamma_s = rho_w * g * gamma_s_over_rho_g;

    K_h_f   = rho_w * g * A_w_f;              % hydrostatic stiffness
    K_h_s   = rho_w * g * A_w_s;
    
    H       = Hs / sqrt(2);                   % equivalent regular wave height

    F_f_mag = gamma_f .* H / 2;               % excitation force of wave
    F_s_mag = gamma_s .* H / 2;

    F_f_phase = gamma_phase_f;                % excitation phase
    F_s_phase = gamma_phase_s;

    A_drag_s = pi/4 * D_d^2;

                                              % values used for drag
    drag_const_f = 4/(3*pi) * rho_w * A_w_f * C_d_float;
    drag_const_s = 4/(3*pi) * rho_w * A_drag_s * C_d_spar;
    mag_v0_f = H/2 * g .* k ./ w .* exp(-k * T_f/2); % finite depth velocity of incident wave
    mag_v0_s = H/2 * g .* k ./ w .* exp(-k * T_s/2);
end

function k = dispersion(w,h,g)
    h_lambda_deep = 0.4; % h/lambda > threshold for deep
    h_lambda_shallow = 0.05; % h/lambda < threshold for shallow
    idx_deep = w > sqrt(h_lambda_deep * 2*pi * g / h); 
    idx_shallow = w < h_lambda_shallow * 2*pi * sqrt(g/h);
    idx_mid = ~idx_deep & ~idx_shallow;

    k = zeros(size(w));
    k(idx_deep) = w(idx_deep).^2 / g;
    k(idx_shallow) = w(idx_shallow) / sqrt(g*h);

    if any(idx_mid,'all')
        err = 1;
        iters = 0;
        k_deep_guess = w(idx_mid).^2 / g;
        k_guess = k_deep_guess;
        while err > 0.2 && iters < 50
            k_new = w(idx_mid).^2 / g ./ tanh(k_guess*h);
            err = max(abs(k_new - k_guess)./k_deep_guess,[],'all');
            k_guess = k_new;
            iters = iters + 1;
        end
        k(idx_mid) = k_new;
    end
end
