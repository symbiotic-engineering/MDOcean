function b = var_bounds(p)

b.D_sft_min = 0;
b.D_sft_max = 50;
b.D_sft_nom = 20;

b.D_i_ratio_min = 0;
b.D_i_ratio_max = 1;
b.D_i_ratio_nom = .3;

b.D_or_ratio_min = 1;
b.D_or_ratio_max = 5;
b.D_or_ratio_nom = 1.5;

b.M_min = 1;
b.M_max = length(p.sigma_y);
b.M_nom = 1;

b.F_max_min = 0;
b.F_max_max = 100;
b.F_max_nom = 1;

b.D_int_min = 1;
b.D_int_max = 100;
b.D_int_nom = 10;

b.w_n_min = 2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
b.w_n_max = 2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
b.w_n_nom = 2*pi/8;

end