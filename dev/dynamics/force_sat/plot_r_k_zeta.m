% plot r_k * zeta as a function of omega for optimal control
% for the nominal design and LCOE-optimal design
clear

p = parameters();
b = var_bounds();

nominalX = [b.X_noms; 1];
optimalX = gradient_optim(b.X_start_struct,p,b,1);
Xs = [nominalX, optimalX];

figure
hold on

linespec = {'','--'};
for i = 1:2
    Xi = Xs(:,i);
    [D_f, D_s_ratio, h_f_ratio, T_s_ratio] = deal(Xi(1),Xi(2),Xi(3),Xi(4));
    h_f = b.h_f_ratio_nom * D_f;
    T_f = p.T_f_over_h_f * h_f;
    D_s = D_s_ratio * D_f;
    T_s = p.T_s_over_D_s * D_s;
    h_s = T_s / b.T_s_ratio_nom;
    h_d = p.h_d_over_D_s * D_s;
    D_d = p.D_d_over_D_s * D_s;

    w = linspace(0.1,2);
    T = 2*pi./w;
    Hs = 1;
    [w,A,B_h,K_h] = dynamics_simple(Hs, T, D_f, T_f, p.rho_w, p.g);
    [~, ~, m_float] = geometry(D_s, D_f, T_f, h_f, h_s, ...
                                                p.t_ft, p.t_fr, p.t_fc, p.t_fb, p.t_sr, ...
                                                p.t_dt, D_d, p.D_dt, p.theta_dt, T_s, h_d, ...
                                                b.M_nom, p.rho_m, p.rho_w, p.m_scale);
    m = m_float + A;

    r_k_zeta = w .* B_h ./ (m*w.^2 - K_h);


    plot(w,abs(r_k_zeta),linespec{i})
end

xlabel('\omega')
ylabel('|r_k \zeta|')
title('Simulated Hydro Coeffs')
ylim([-.5 3])
grid on
legend('Nominal','Optimal LCOE')
improvePlot
