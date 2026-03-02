n = 5; % number of points per dimension

zero_to_one   = linspace(0.01,0.99,n);
zero_to_three = linspace(0.01,3,   n);

m0h   = zero_to_three;
a1_a2 = zero_to_one;
m0a2  = zero_to_three;
d1_h  = zero_to_one;
d2_d1 = zero_to_one;
a3_a2 = zero_to_three;

[M0H, A1_A2, M0A2, ...
 D1_H, D2_D1, A3_A2] = ndgrid(m0h, a1_a2, m0a2, d1_h, d2_d1, a3_a2);

A2_H = M0A2  ./ M0H;
A1_H = A1_A2 .* A2_H;
D2_H = D2_D1 .* D1_H;
A3_H = A3_A2 .* A2_H;

p = parameters();

A1 = A1_H * p.h;
A2 = A2_H * p.h;
A3 = A3_H * p.h;
D1 = D1_H * p.h;
D2 = D2_H * p.h;

b = var_bounds();
X = [b.X_noms; 1];

for i=1:length(A1)
    D_s = 2*A1(i);
    D_f = 2*A2(i);
    D_d = 2*A3(i);
    T_s = D1;
    T_f_2 = D2;

    % construct params and design vec
    p_i = p;
    p_i.T_s_over_D_s = T_s / D_s;
    p_i.D_d_over_D_s = D_d / D_s;
    X_i = X;
    X_i(strcmp(b.var_names),'D_s') = D_s;
    X_i(strcmp(b.var_names),'D_f') = D_f;
    X_i(strcmp(b.var_names),'h_s') = T_s + 5;
    X_i(strcmp(b.var_names),'T_f_2') = T_f_2;
    
    % run simulation (drag and force/power sat are turned off inside check_max_CW)
    [hydro_ratio, ~, ...
     ~, ~, ~, ~, ~, val] = check_max_CW('', p_i, X_i, false);
end

% NMK = 10;
% H = 1;
% 
% g_nondim = g; % fixme so mult works out
% w_nondim = w; % fixme so mult works out
% 
% % mult = M0H .* (1 - w2_over_gk.^2) + w2_over_gk;
% C_d_float = 0; C_d_spar = 0;
% use_MEEM = true;
% spar_excitation_coeffs = get_spar_exc(g);
% hydro = []; % unused (only for use_MEEM = false)
% 
% wave_amp = 1;
% wave_height = 2*wave_amp;
% Hs = sqrt(2)*wave_height;
% 
% [m_f,B_h_f,K_h_f,F_f_mag,F_f_phase,...
%  m_s,B_h_s,K_h_s,F_s_mag,F_s_phase,...
%  m_c,B_h_c,drag_const_f,drag_const_s,...
%  mag_v0_f,mag_v0_s,w,k,A_f_over_rho, ...
%  A_s_over_rho, A_c_over_rho, ...
%  B_f_over_rho_w, B_s_over_rho_w, ...
%  B_c_over_rho_w, gamma_f_over_rho_g, ...
%  gamma_s_over_rho_g, gamma_phase_f, ...
%                      gamma_phase_s] = get_dynamic_coeffs(Hs, T, D_f, T_f, ...
%                                                         D_s, D_d, T_s, h, ...
%                                                         m_float, m_spar, ...
%                                                         spar_excitation_coeffs,...
%                                                         C_d_float, C_d_spar, ...
%                                                         rho_w, g, use_MEEM, ...
%                                                         harmonics, hydro);
% [mag_U,phase_U,...
%  real_P,reactive_P,...
%  mag_X_u,phase_X_u,...
%  mag_X_f,phase_X_f,...
%  mag_X_s,phase_X_s,...
%  B_p,K_p,P_sat_ratio] = get_response_drag(w,m_f,m_s,m_c,...
%                                 B_h_f,B_h_s,B_h_c,...
%                                 K_h_f,K_h_s,...
%                                 F_f_mag,F_f_phase,F_s_mag,F_s_phase,F_max,...
%                                 drag_const_f,drag_const_s,mag_v0_f,mag_v0_s,...
%                                 X_max,control_type,multibody,...
%                                 X_tol,phase_X_tol,max_drag_iters,...
%                                 drag_convergence_plot_on);
