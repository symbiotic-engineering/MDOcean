%% Updated Equations
syms D_s h_s h_f D_f rho_w g depth A_f_l T_s t_fc T_f real positive
syms A_d_l P_hydrostatic A_f t_sr V_rims_gussets h_d A_f_c real positive
syms A_vc_c A_d_c A_vc_l D_d CG_vc num_gussets_loaded_lateral real positive
syms sigma_zz1_temp D_s_temp V_vc_m t_fr num_gussets real positive
syms D_vc_i real positive
syms principal_term [3,1] real
syms shear_term [3,1] real
syms s_vm [3,1] real positive
syms sigma_rr [3,1] real positive
syms sigma_tt [3,1] real positive
syms sigma_zz [3,1] real positive
syms sigma_tz [3,1] real positive
syms sigma_rt [3,1] real positive
syms sigma_zr [3,1] real positive
syms F_surge [3,1] real positive
syms F_heave [3,1] real positive

% set up the von Mises stress equations
principal_term(1) = 1/2 * ((sigma_rr(1) - sigma_tt(1)).^2 + (sigma_tt(1) - sigma_zz(1)).^2 ...
    + (sigma_zz(1) - sigma_rr(1)).^2);
principal_term(2) = 1/2 * ((sigma_rr(2) - sigma_tt(2)).^2 + (sigma_tt(2) - sigma_zz(2)).^2 ...
    + (sigma_zz(2) - sigma_rr(2)).^2);
principal_term(3) = 1/2 * ((sigma_rr(3) - sigma_tt(3)).^2 + (sigma_tt(3) - sigma_zz(3)).^2 ...
    + (sigma_zz(3) - sigma_rr(3)).^2);
shear_term(1) = 3 * (sigma_rt(1).^2 + sigma_tz(1).^2 + sigma_zr(1).^2);
shear_term(2) = 3 * (sigma_rt(2).^2 + sigma_tz(2).^2 + sigma_zr(2).^2);
shear_term(3) = 3 * (sigma_rt(3).^2 + sigma_tz(3).^2 + sigma_zr(3).^2);

eqn1 = s_vm(1) == sqrt(principal_term(1) + shear_term(1));
eqn2 = s_vm(2) == sqrt(principal_term(2) + shear_term(2));
eqn3 = s_vm(3) == sqrt(principal_term(3) + shear_term(3));

% substitutions are commented
% A_lat_sub = [A_f_l A_vc_l A_d_l] -> sigma_surge = F_surge ./ A_lat_sub ->
% sigma_rr = P_hydrostatic + sigma_surge
eqn4 = sigma_rr(1) == P_hydrostatic + (F_surge(1) / A_f_l);
% r_over_t = [0, D_s/(2*t_sr), 0] -> sigma_tt = P_hydrostatic .* r_over_t
eqn5 = sigma_tt(1) == P_hydrostatic;
% V_rims_gussets = A_f_c * h_f -> A_c = [A_f_c, A_vc_c, A_d_c] ->
% sigma_zz = F_heave ./ A_c
eqn6 = sigma_zz1_temp == F_heave(1) * h_f / V_rims_gussets;

%eqn1 = subs(eqn1, sigma_zz1, F_heave(1) * h_f / V_rims_gussets)
eqn1 = subs(eqn1, sigma_zz1, F_heave(1) ./ A_f_c);
eqn1 = subs(eqn1, A_f_c, V_rims_gussets/h_f);
% there is a plus and minus in h_f; verify which one is correct
h_f = solve(eqn1, h_f);
h_f = h_f(1);
% A_lat_sub = [A_f_l A_vc_l A_d_l] -> sigma_surge = F_surge ./ A_lat_sub ->
% sigma_rt = sigma_surge
eqn7 = sigma_rt(1) == F_surge(1) / A_f_l;
h_f = subs(h_f, sigma_tz(1), 0);
%eqn9 =  sigma_zr(1) == 0;
h_f = subs(h_f, sigma_zr(1), 0);
h_f = subs(h_f, sigma_rr(1), P_hydrostatic + (F_surge(1) / A_f_l));
h_f = subs(h_f, sigma_tt(1), P_hydrostatic);
h_f = subs(h_f, sigma_rt(1), F_surge(1) / A_f_l);
h_f = subs(h_f, A_f_l, num_gussets_loaded_lateral * t_fc * T_f)

eqn10 = sigma_rr(2) == P_hydrostatic + (F_surge(2) / (1/2 * pi * D_s * T_s));
eqn11 = sigma_tt(2) == P_hydrostatic + (sqrt(D_f^2-(4*(pi/4 * (D_f^2 - D_s^2)) / pi)))/(2 * t_sr);
eqn12 = sigma_zz(2) == F_heave(2) / A_vc_c;
eqn13 = sigma_rt(2) == F_surge(2) / A_vc_l;
eqn14 = sigma_tz(2) == 0;
eqn15 =  sigma_zr(2) == 0;

% solve for D_s
eqn2 = subs(eqn2, sigma_rt(2), F_surge(2) / (1/2 * pi * D_s * T_s));
eqn2 = subs(eqn2, sigma_tt(2), P_hydrostatic + (D_s/(2*t_sr)));
D_s = solve(eqn2, D_s);
% D_s has a plus or minus; verify which one is correct
D_s = D_s(1);
D_s = subs(D_s, sigma_zz(2), F_heave(2) / V_vc_m / h_s);
D_s = subs(D_s, V_vc_m, A_vc_c * h_s);
D_s = subs(D_s, A_vc_c, pi/4 * (D_s^2 - D_vc_i^2));
D_s = subs(D_s, sigma_rt(2), F_surge(2) / (1/2 * pi * D_s * T_s));
D_s = subs(D_s, sigma_tz(2), 0);
D_s = subs(D_s, sigma_rr(2), P_hydrostatic + (F_surge(2) ./ A_vc_l));
D_s = subs(D_s, sigma_zr(2), 0)

% A_lat_sub = [A_f_l A_vc_l A_d_l] -> sigma_surge = F_surge ./ A_lat_sub ->
% sigma_rr = P_hydrostatic + sigma_surge
eqn16 = sigma_rr(3) == P_hydrostatic + (F_surge(3) / A_d_l);
eqn17 = sigma_tt(3) == P_hydrostatic;
eqn18 = sigma_zz(3) == F_heave(3) / A_d_c;
eqn19 = sigma_rt(3) == F_surge(3) / (1/2 * pi * D_d * (CG_vc - (h_s/2)));
eqn20 = sigma_tz(3) == 0;
eqn21 = sigma_zr(3) == 0;

eqn3 = subs(eqn3, sigma_rt(3), F_surge(3) / (1/2 * pi * D_d * (CG_vc - (h_s/2))));
D_d = solve(eqn3, D_d);
% there is a plus and minus; verify which one is correct
D_d = D_d(1);
D_d = subs(D_d, sigma_tt(3), P_hydrostatic);
D_d = subs(D_d, sigma_zz(3), F_heave(3) / (pi/4 * (D_d^2 - D_s^2)));
D_d = subs(D_d, sigma_rt(3), F_surge(3) / (1/2 * pi * D_d * (CG_vc - (h_s/2))));
D_d = subs(D_d, sigma_tz(3),0);
D_d = subs(D_d, sigma_zr(3),0);
D_d = subs(D_d, sigma_rr(3), P_hydrostatic + (F_surge(3) / (1/2 * pi * D_d * h_d)));
D_d = subs(D_d, CG_vc, h_d + h_s/2)