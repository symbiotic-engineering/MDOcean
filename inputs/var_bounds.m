function b = var_bounds(p)

b.D_f_min = 6;    % because D_s_nom = 6, a consequence of the natural frequency constraint
b.D_f_max = 50;
b.D_f_nom = 20;

b.D_s_ratio_min = 0;
b.D_s_ratio_max = 0.99;
b.D_s_ratio_nom = 6/20;

b.h_f_ratio_min = .1;
b.h_f_ratio_max = 10;
b.h_f_ratio_nom = 4/20;

b.T_s_ratio_min = 0;
b.T_s_ratio_max = 0.99;
b.T_s_ratio_nom = 35/44;

b.F_max_min = 1e-3;
b.F_max_max = 100;
b.F_max_nom = 1;

b.D_int_min = .1;
b.D_int_max = 100;
b.D_int_nom = 10;

b.w_n_min = .01;%2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
b.w_n_max = 41;%2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
b.w_n_nom = 2*pi/8;

b.M_min = 1;
b.M_max = length(p.sigma_y);
b.M_nom = 1;

b.X_mins = [b.D_f_min b.D_s_ratio_min b.h_f_ratio_min b.T_s_ratio_min b.F_max_min b.D_int_min b.w_n_min]';
b.X_maxs = [b.D_f_max b.D_s_ratio_max b.h_f_ratio_max b.T_s_ratio_max b.F_max_max b.D_int_max b.w_n_max]';
b.X_noms = [b.D_f_nom b.D_s_ratio_nom b.h_f_ratio_nom b.T_s_ratio_nom b.F_max_nom b.D_int_nom b.w_n_nom]';

end