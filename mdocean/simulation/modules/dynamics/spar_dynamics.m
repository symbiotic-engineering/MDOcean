
function [A_s_over_rho,gamma_s_over_rho_g,...
    B_s_over_rho_w,gamma_s_phase,A_c_over_rho,B_c_over_rho_w] = spar_dynamics(Ds_over_Dd, D_d, T_s, ...
                                                    spar_coeffs, k, mult)

    % added mass
    r = Ds_over_Dd;
    root_r = sqrt(1-r^2);
    ratio_term = 1/3 - 1/4 * r^2 * root_r - 1/12 * (1 - root_r)^2 * (2 + root_r);
    A_s_over_rho = D_d^3 * ratio_term; % added mass

    % excitation - interpolate WAMIT results
    depth_multiplier = exp(-k * (T_s - spar_coeffs.T_s));
    gamma_s_over_rho_g = interp1(spar_coeffs.k * D_d, spar_coeffs.gamma_over_rho_g, k * D_d) .* depth_multiplier;
    gamma_s_phase      = interp1(spar_coeffs.k * D_d, spar_coeffs.gamma_phase,      k * D_d);
    A_c_over_rho       = interp1(spar_coeffs.k * D_d, spar_coeffs.A_c_over_rho,     k * D_d);
    B_c_over_rho_w     = interp1(spar_coeffs.k * D_d, spar_coeffs.B_c_over_rho_w,   k * D_d);

    % radiation damping
    B_s_over_rho_w = k/2 .* gamma_s_over_rho_g.^2 ./ mult; % finite depth haskind relation

%     % frequency parameter
%     f = w / (2 * pi); % frequency Hz
%     nu = mu / rho_w; % kinematic viscosity of water
%     beta = D_d^2 * f / nu; % nondimensional frequency parameter
%
%     x_s_error = ones(size(w));
%     x_s_guess = ones(size(w));
%     tol = 0.01;
%     num_iter = 0;
%     while any( abs(x_s_error) > tol, 'all')
%         [x_s, x_s_error] = calculate_x_s(x_s_guess, D_d, C_d, K_s, K_p, m_s, ...
%                                         B_p, B_33_rad_spar, w, x_f, F_s, mu, beta);
%         x_s_guess = x_s;
%         num_iter = num_iter + 1;
%     end
%     fprintf('Took %d iterations to converge \n',num_iter)

    %x_s_x_f_ratio = x_s ./ x_f
    %Kp_over_Ks = K_p ./ K_s;
end

% todo: test spar equations, make figure of wamit results, latex equations

function [x_s, x_s_error] = calculate_x_s(x_s_guess, D_d, C_d, K_s, K_p, m_s, B_p, B_33_rad_spar, w, x_f, F_s, mu, beta)

    % drag damping
    KC = 2 * pi * x_s_guess / D_d; % fixme check if this is D_d or D
    B_33_drag_spar = 1/3 * mu * beta * D_d .* KC * C_d;
    B_s = B_33_drag_spar + B_33_rad_spar; % total damping from drag and radiation

    % transfer functions
    s = 1i * w;
    x_s_over_F_s = 1./ (K_p + K_s + (B_p+B_s) .* s + m_s .* s.^2);
    x_s_over_x_f = (K_p + B_p .* s) .* x_s_over_F_s;

    x_s_complex = x_s_over_F_s .* F_s + x_s_over_x_f .* x_f; % fixme I think x_f should be complex since it may have different phase than F_s
    x_s = abs(x_s_complex);

    x_s_error = x_s_guess - x_s;

    %zeta_with_coupling = (B_p+B_s)./(2*sqrt((K_p+K_s).*m_s))

    %if any( ~isreal(zeta_with_coupling),'all')
    %    disp('ohno')
    %end
    %zeta_no_coupling = B_s./(2*sqrt(K_s.*m_s))

end
