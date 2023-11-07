% %function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling,D_s,D_f,h_s,t_ft,T_f,h_f_1,t_fb,T_s,h_d]=v2_von_Mises_Symbolic(F_heave,F_surge,M,h_s,T_s,rho_w,g,sigma_y,A_c,A_lat_sub,r_over_t,I,E)
% 
% syms FOS t_fc t_fb CB_vc real positive
% syms D_f t_fr h_d h_s F_heave g rho_w F_surge real positive
% syms A_c A_f_c A_f_l t_sr D_s m_f_b m_f_tot m_f_m real positive
% syms T_s T_f r_over_t V_f_b  V_f_pct mass real positive
% syms sigma_rr sigma_tt sigma_zz sigma_rt sigma_tz sigma_zr real positive
% syms A_f V_vc_m V_f_d V_f_tot P_hydrostatic sigma_surge real positive
% syms h_f V_top_plate t_ft V_bot_plate V_vc_d real positive
% syms A_lat_sub [1,3] real positive
% syms I [1,3] real positive
% syms sigma_y [1,3] real positive
% syms s_vm [3,3] real positive
% syms mass [3,1] real positive
% syms V_d [3,1] real positive
% 
% 
% E=[1 1 1];
% M=1;
% %%%%%%%% update these
% num_gussets = 24;
% num_gussets_loaded_lateral = 2;
% %depth = [0 T_s T_s];
% %P_hydrostatic = rho_w * g * depth;
% %sigma_surge = F_surge ./ A_lat_sub;
% 
% % get the equation for s_vm
% principal_term = 1/2 * ((sigma_rr - sigma_tt).^2 + (sigma_tt - sigma_zz).^2 ...
%     + (sigma_zz - sigma_rr).^2);
% shear_term = 3 * (sigma_rt.^2 + sigma_tz.^2 + sigma_zr.^2);
% eqn = s_vm == sqrt(principal_term + shear_term);
% for i = 1:3
%     for j = 1:3
%         sigma_rr_temp = solve(eqn(i,j), sigma_rr);
%         sigma_rr_temp2(i,j) = sigma_rr_temp(1);
%     end
% end
% sigma_rr = sigma_rr_temp2;
% 
% % substitute in stresses
% sigma_rr = subs(sigma_rr, sigma_tt, P_hydrostatic .* r_over_t);
% sigma_rr = subs(sigma_rr, sigma_zz, F_heave ./ A_c);
% sigma_rr = subs(sigma_rr, sigma_rt, sigma_surge);
% sigma_rr = subs(sigma_rr, sigma_tz, 0);
% sigma_rr = subs(sigma_rr, sigma_zr, 0)
% eqn_sigma_rr = sigma_rr == P_hydrostatic + sigma_surge;
% 
% % find intermediate variables
% for i = 1:3
%     for j = 1:3
%         r_over_t_temp = solve(eqn_sigma_rr(i,j), r_over_t);
%         r_over_t_temp2(i,j) = r_over_t_temp(1);
%     end
% end
% r_over_t_1 = r_over_t_temp2;
% 
% % D_s
% % in terms of t_sr, F_heave, A_c, sigma_surge, P_hydrostatic, s_vm
% for i = 1:3
%     for j = 1:3
%         D_s_temp = solve(r_over_t_1(i,j) == D_s/(2*t_sr),D_s);
%         D_s_temp2(i,j) = D_s_temp(1);
%     end
% end
% D_s = D_s_temp2;
% 
% % D_f
% % in terms of t_sr, F_heave, A_c, sigma_surge, P_hydrostatic, s_vm, A_f
% for i = 1:3
%     for j = 1:3
%         D_f_temp = solve(A_f == pi/4 * (D_f^2 - (D_s(i,j))^2), D_f);
%         D_f_temp2(i,j) = D_f_temp;
%     end
% end
% D_f = D_f_temp2
% 
% % h_s
% % in terms of t_sr, F_heave, A_c, sigma_surge, P_hydrostatic, s_vm, V_vc_m
% D_vc_i = D_s - 2 * t_sr;
% A_vc_c = pi/4 * (D_s.^2 - D_vc_i.^2);
% for i = 1:3
%     for j = 1:3
%         h_s_temp = solve(V_vc_m == A_vc_c(i,j) * h_s, h_s);
%         h_s_temp2(i,j) = h_s_temp;
%     end
% end
% h_s = h_s_temp2;
% 
% % h_f
% % in terms of A_f, rho_w, mass, V_d, V_f_pct
% eqn_h_f = V_f_tot == A_f * h_f;
% h_f = solve(eqn_h_f, h_f);
% h_f = subs(h_f, V_f_tot, V_f_b / V_f_pct);
% h_f = subs(h_f, V_f_b, m_f_b / rho_w);
% h_f = subs(h_f, m_f_b, m_f_tot - m_f_m);
% h_f = subs(h_f, m_f_tot, V_f_d * rho_w);
% h_f = subs(h_f, m_f_m, mass(1));
% h_f = subs(h_f, V_f_d, V_d(1));
% 
% % t_ft
% % in terms of t_sr, F_heave, A_c, sigma_surge, P_hydrostatic, s_vm, A_f,
% % V_top_plate
% for i = 1:3
%     for j = 1:3
%         t_ft_temp = solve(V_top_plate == pi * (D_f(i,j)./2)^2 * t_ft, t_ft);
%         t_ft_temp2(i,j) = t_ft_temp;
%     end
% end
% t_ft = t_ft_temp2
% 
% % t_fb
% % in terms of 
% for i = 1:3
%     for j = 1:3
%         t_fb_temp = solve(V_bot_plate == pi * (D_f(i,j)./2)^2 * t_fb, t_fb);
%         t_fb_temp2(i,j) = t_fb_temp;
%     end
% end
% t_fb = t_fb_temp2;
% 
% % % T_s
% % eqn_T_s = V_vc_d == pi/4 * D_s1^2 * T_s;
% % T_s = solve(eqn_T_s, T_s);
% for i = 1:3
%     for j = 1:3
%         T_s_temp = solve(V_vc_d == pi/4 * D_s(i,j)^2 * T_s, T_s);
%         T_s_temp2(i,j) = T_s_temp;
%     end
% end
% T_s = T_s_temp2;
% 
% % % T_f
% % eqn_t_fr = A_f_c == pi * (D_f1 + D_s1) * t_fr + num_gussets * t_fc * (D_f1 - D_s1)/2;
% % t_fr = solve(eqn_t_fr, t_fr);
% % eqn_t_fc = A_f_l == num_gussets_loaded_lateral * t_fc * T_f;
% % t_fc_new = solve(eqn_t_fc, t_fc);
% % t_fr = subs(t_fr, t_fc, t_fc_new);
% % eqn_T_f = V_f_d == A_f * T_f;
% % T_f_new = solve(eqn_T_f, T_f);
% % t_fr = subs(t_fr, T_f, T_f_new);
% 
% % % h_d
% % eqn_h_d = CB_vc == h_d + T_s/2;
% % h_d = solve(eqn_h_d, h_d);
% % 
% % K = 2; % fixed-free - top is fixed by float angular stiffness, bottom is free
% % L = h_s;
% % F_buckling = pi^2 * E(M) * I(2) / (K*L)^2;
% % 
% % % Factor of Safety (FOS) Calculations
% % FOS_yield = sigma_y(M) ./ s_vm;
% % FOS1Y = FOS_yield(1);
% % FOS2Y = FOS_yield(2);
% % FOS3Y = FOS_yield(3);
% % FOS_buckling = F_buckling ./ F_heave;
% 
% filename = 'simulation\modules\dynamics\von_Mises_Symbolic';
% %matlabFunction(D_f, D_s, h_f, h_s,'File',filename,'Vars',[F_heave,F_surge,M,h_s,T_s,rho_w,g,sigma_y,A_c,A_lat_sub,r_over_t,I,E,mass]);
% %matlabFunction([D_f, D_s, h_f, h_s],"File",filename)
% % matlabFunction(D_f, D_s, h_f, h_s,'File',filename)


%% Updated Equations
syms D_s h_s h_f D_f rho_w g depth A_f_l T_s real positive
syms A_d_l P_hydrostatic A_f t_sr V_rims_gussets real positive
syms A_vc_c A_d_c A_vc_l D_d CG_vc real positive
syms principal_term [3,1] real
syms shear_term [3,1] real
syms s_vm [3,1] real
syms sigma_rr [3,1] real positive
syms sigma_tt [3,1] real positive
syms sigma_zz [3,1] real positive
syms sigma_tz [3,1] real positive
syms sigma_rt [3,1] real positive
syms sigma_zr [3,1] real positive
syms F_surge [3,1] real positive
syms F_heave [3,1] real positive

principal_term(1) = 1/2 * ((sigma_rr(1) - sigma_tt(1)).^2 + (sigma_tt(1) - sigma_zz(1)).^2 ...
    + (sigma_zz(1) - sigma_rr(1)).^2);
principal_term(2) = 1/2 * ((sigma_rr(2) - sigma_tt(2)).^2 + (sigma_tt(2) - sigma_zz(2)).^2 ...
    + (sigma_zz(2) - sigma_rr(2)).^2);
principal_term(3) = 1/2 * ((sigma_rr(3) - sigma_tt(3)).^2 + (sigma_tt(3) - sigma_zz(3)).^2 ...
    + (sigma_zz(3) - sigma_rr(3)).^2);
shear_term(1) = 3 * (sigma_rt(1).^2 + sigma_tz(1).^2 + sigma_zr(1).^2);
shear_term(2) = 3 * (sigma_rt(2).^2 + sigma_tz(2).^2 + sigma_zr(2).^2);
shear_term(3) = 3 * (sigma_rt(3).^2 + sigma_tz(3).^2 + sigma_zr(3).^2);

% Make it so body 1 has the first entry of any vector, body 2 have the
% second entry, etc.
% Solve using all three bodies, not one
% Solve for D_s, h_s, and h_f in terms of D_f

% Remove equations from geometry that have alreay been used
% In simulation module, remove D_s over ... (There are three of them?)
% Run using gradientoptim.m

eqn1 = s_vm(1) == sqrt(principal_term(1) + shear_term(1));
eqn2 = s_vm(2) == sqrt(principal_term(2) + shear_term(2));
eqn3 = s_vm(3) == sqrt(principal_term(3) + shear_term(3));
% A_lat_sub = [A_f_l A_vc_l A_d_l] -> sigma_surge = F_surge ./ A_lat_sub ->
% sigma_rr = P_hydrostatic + sigma_surge
eqn4 = sigma_rr(1) == P_hydrostatic + (F_surge(1) / A_f_l);
% r_over_t = [0, D_s/(2*t_sr), 0] -> sigma_tt = P_hydrostatic .* r_over_t
eqn5 = sigma_tt(1) == P_hydrostatic;
% V_rims_gussets = A_f_c * h_f -> A_c = [A_f_c, A_vc_c, A_d_c] ->
% sigma_zz = F_heave ./ A_c
eqn6 = sigma_zz(1) == F_heave(1) / (V_rims_gussets/h_f);
% A_lat_sub = [A_f_l A_vc_l A_d_l] -> sigma_surge = F_surge ./ A_lat_sub ->
% sigma_rt = sigma_surge
eqn7 = sigma_rt(1) == F_surge(1) / A_f_l;
eqn8 = sigma_tz(1) == 0;
eqn9 =  sigma_zr(1) == 0;
eqn10 = sigma_rr(2) == P_hydrostatic + (F_surge(2) / (1/2 * pi * D_s * T_s));
eqn11 = sigma_tt(2) == P_hydrostatic + (sqrt(D_f^2-(4*A_f / pi)))/(2 * t_sr);
eqn12 = sigma_zz(2) == F_heave(2) / A_vc_c;
eqn13 = sigma_rt(2) == F_surge(2) / A_vc_l;
eqn14 = sigma_tz(2) == 0;
eqn15 =  sigma_zr(2) == 0;
% A_lat_sub = [A_f_l A_vc_l A_d_l] -> sigma_surge = F_surge ./ A_lat_sub ->
% sigma_rr = P_hydrostatic + sigma_surge
eqn16 = sigma_rr(3) == P_hydrostatic + (F_surge(3) / A_d_l);
eqn17 = sigma_tt(3) == P_hydrostatic;
eqn18 = sigma_zz(3) == F_heave(3) / A_d_c;
eqn19 = sigma_rt(3) == F_surge(3) / (1/2 * pi * D_d * (CG_vc - (h_s/2)));
eqn20 = sigma_tz(3) == 0;
eqn21 = sigma_zr(3) == 0;

eqns = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10, ...
    eqn11, eqn12, eqn13, eqn14, eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21];

s = symvar(eqns)

unknowns = [sigma_rr(1), sigma_rt(1), sigma_tt(1), sigma_tz(1), sigma_zr(1), ...
    sigma_zz(1), sigma_rr(2), sigma_rt(2), sigma_tt(2), sigma_tz(2), ...
    sigma_zr(2), sigma_zz(2), sigma_rr(3), sigma_rt(3), sigma_tt(3), ...
    sigma_tz(3), sigma_zr(3), sigma_zz(3), D_s, h_s, h_f];

[sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9, sol10, sol11, sol12, ...
    sol13, sol14, sol15, sol16, sol17, sol18, sol19, sol20, sol21, params, conds] ...
    = solve(eqns,unknowns,'ReturnConditions',true,'IgnoreAnalyticConstraints',true)