function b = var_bounds(p)

b.D_sft_min = 0;
b.D_sft_max = 50;
b.D_sft_nom = 20;

b.D_i_ratio_min = 0;
b.D_i_ratio_max = 1;
b.D_i_ratio_nom = .3;

b.D_or_min = 0;
b.D_or_max = 50;
b.D_or_nom = 30;

b.M_min = 1;
b.M_max = length(p.sigma_y);
b.M_nom = 1;

b.N_WEC_min = 1;
b.N_WEC_max = 100;
b.N_WEC_nom = 10;

b.D_int_min = 1;
b.D_int_max = 100;
b.D_int_nom = 10;

b.w_n_min = 2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
b.w_n_max = 2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
b.w_n_nom = 2*pi/8;

end