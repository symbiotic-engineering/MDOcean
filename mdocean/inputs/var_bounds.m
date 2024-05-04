function b = var_bounds()

b.D_f_min = 6;    % because D_s_nom = 6, a consequence of the spar natural frequency constraint
b.D_f_max = 40;
b.D_f_nom = 20;
b.D_f_start = 20;

b.D_s_over_D_f_min = 0.01;
b.D_s_over_D_f_max = 0.99;
b.D_s_over_D_f_nom = 6/20;
b.D_s_over_D_f_start = 6/20;

b.h_f_over_D_f_min = .1;
b.h_f_over_D_f_max = 10;
b.h_f_over_D_f_nom = 4/20;
b.h_f_over_D_f_start = 4/20;

b.T_s_over_h_s_min = 0.01;
b.T_s_over_h_s_max = 0.99;
b.T_s_over_h_s_nom = 35/44;
b.T_s_over_h_s_start = 35/44;

b.F_max_min = 0.01;
b.F_max_max = 100;
b.F_max_nom = 5;
b.F_max_start = 5;

b.B_p_min = .1;
b.B_p_max = 50;
b.B_p_nom = 10;
b.B_p_start = 0.5;

b.w_n_min = .01;%2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
b.w_n_max = 40;%2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
b.w_n_nom = 0.8;
b.w_n_start = 0.8;

b.M_min = 1;
b.M_max = 3;
b.M_nom = 1;
b.M_start = 1;

b.t_ft_min = 0.3 ;
b.t_ft_max = 0.5 ;
b.t_ft_nom = 0.5;
b.t_ft_start = 0.4;

b.t_fr_min = 0.1; 
b.t_fr_max = 0.6;
b.t_fr_nom = 0.44;
b.t_fr_start = 0.2;

b.t_fc_min = 0.1 ; 
b.t_fc_max = 0.6; 
b.t_fc_nom = 0.44; 
b.t_fc_start = 0.2; 

b.t_fb_min = 0.1;
b.t_fb_max = 0.9 ;
b.t_fb_nom = 0.56;
b.t_fb_start = 0.2;

b.t_sr_min = 0.2;
b.t_sr_max = 2.0;
b.t_sr_nom = 1.00;
b.t_sr_start = 0.3;

b.t_dt_min = 0.2;
b.t_dt_max = 2.0;
b.t_dt_nom = 1.00 ;
b.t_dt_start = 0.3;

b.power_max_min = 3;
b.power_max_max = 40;
b.power_max_nom = 20;
b.power_max_start = 5;

b.X_mins = [b.D_f_min b.D_s_over_D_f_min b.h_f_over_D_f_min b.T_s_over_h_s_min b.F_max_min b.B_p_min b.w_n_min b.M_min b.t_ft_min b.t_fr_min b.t_fc_min b.t_fb_min b.t_sr_min b.t_dt_min b.power_max_min]';
b.X_maxs = [b.D_f_max b.D_s_over_D_f_max b.h_f_over_D_f_max b.T_s_over_h_s_max b.F_max_max b.B_p_max b.w_n_max b.M_max b.t_ft_max b.t_fr_max b.t_fc_max b.t_fb_max b.t_sr_max b.t_dt_max b.power_max_min]';
b.X_noms = [b.D_f_nom b.D_s_over_D_f_nom b.h_f_over_D_f_nom b.T_s_over_h_s_nom b.F_max_nom b.B_p_nom b.w_n_nom b.M_nom b.t_ft_nom b.t_fr_nom b.t_fc_nom b.t_fb_nom b.t_sr_nom b.t_dt_nom b.power_max_min]';
b.X_starts = [b.D_f_start b.D_s_over_D_f_start b.h_f_over_D_f_start b.T_s_over_h_s_start b.F_max_start b.B_p_start b.w_n_start b.M_start b.t_ft_start b.t_fr_start b.t_fc_start b.t_fb_start b.t_sr_start b.t_dt_start b.power_max_min]';

b.X_start_struct = struct('D_f',b.D_f_start,'D_s_over_D_f',b.D_s_over_D_f_start,...
        'h_f_over_D_f',b.h_f_over_D_f_start,'T_s_over_h_s',b.T_s_over_h_s_start,...
        'F_max',b.F_max_start,'B_p',b.B_p_start,'w_n',b.w_n_start,'M',b.M_start,...
        't_ft',b.t_ft_start, 't_fr',b.t_fr_start, 't_fc', b.t_fc_start, 't_fb', b.t_fb_start,...
        't_sr', b.t_sr_start, 't_dt', b.t_dt_start, 'power_max', b.power_max_max);

b.var_names = {'D_f','D_s_over_D_f','h_f_over_D_f','T_s_over_h_s','F_max','B_p','w_n','M','t_ft','t_fr','t_fc','t_fb','t_sr','t_dt','power_max'};
b.constraint_names = {'float_too_heavy','float_too_light','spar_too_heavy','spar_too_light',...
    'stability','FOS_float_yield','FOS_col_yield','FOS_plate','FOS_col_buckling',...
    'pos_power','spar_damping','spar_height','LCOE_max','irrelevant_max_force'};

% modify nominal control inputs to so power and force matches actual
[F_max_nom, B_p_nom, w_n_nom] = find_nominal_inputs(b, false);
b.F_max_nom = F_max_nom;
b.B_p_nom = B_p_nom;
b.w_n_nom = w_n_nom;

b.X_noms = [b.D_f_nom b.D_s_over_D_f_nom b.h_f_over_D_f_nom b.T_s_over_h_s_nom b.F_max_nom b.B_p_nom b.w_n_nom b.M_nom b.t_ft_nom b.t_fr_nom b.t_fc_nom b.t_fb_nom b.t_sr_nom b.t_dt_nom b.power_max_nom]';

end