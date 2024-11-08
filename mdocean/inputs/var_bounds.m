function b = var_bounds()

b.var_names = {'D_f','D_s_ratio','h_f_ratio','T_s_ratio','F_max','D_int','w_n','M'};
b.var_names_pretty = {'D_f','D_s/D_f','h_f/D_f','T_s/h_s','F_{max}','B_p','\omega_n','M'};

% outer diameter of float (m)	
b.D_f_min = 6;    % because D_s_nom = 6, a consequence of the spar natural frequency constraint
b.D_f_max = 40;
b.D_f_nom = 20;
b.D_f_start = 20;

% inner diameter to outer diameter of float ratio (-)	
b.D_s_ratio_min = 0.01;
b.D_s_ratio_max = 0.99;
b.D_s_ratio_nom = 6/20;
b.D_s_ratio_start = 6/20;

% height of float to outer diameter of float ratio (-)
b.h_f_ratio_min = .1;
b.h_f_ratio_max = 10;
b.h_f_ratio_nom = 5/20;
b.h_f_ratio_start = 5/20;

% draft of spar to height of spar ratio (-)
b.T_s_ratio_min = 0.01;
b.T_s_ratio_max = 0.99;
b.T_s_ratio_nom = 29/37.9;%35/44;
b.T_s_ratio_start = 29/37.9;%35/44;

% maximum powertrain force (MN)
b.F_max_min = 0.01;
b.F_max_max = 100;
b.F_max_nom = 5;
b.F_max_start = 5;

% powertrain damping (MN / (m/s))
b.B_p_min = .1;
b.B_p_max = 50;
b.B_p_nom = 10;
b.B_p_start = 0.5;

% natural frequency (rad/s)
b.w_n_min = .01;%2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
b.w_n_max = 40;%2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
b.w_n_nom = 0.8;
b.w_n_start = 0.8;

% material index (-)
b.M_min = 1;
b.M_max = 3;
b.M_nom = 1;
b.M_start = 1;

b.X_mins = [b.D_f_min b.D_s_ratio_min b.h_f_ratio_min b.T_s_ratio_min b.F_max_min b.B_p_min b.w_n_min]';
b.X_maxs = [b.D_f_max b.D_s_ratio_max b.h_f_ratio_max b.T_s_ratio_max b.F_max_max b.B_p_max b.w_n_max]';
b.X_noms = [b.D_f_nom b.D_s_ratio_nom b.h_f_ratio_nom b.T_s_ratio_nom b.F_max_nom b.B_p_nom b.w_n_nom]';
b.X_starts = [b.D_f_start b.D_s_ratio_start b.h_f_ratio_start b.T_s_ratio_start b.F_max_start b.B_p_start b.w_n_start]';

b.X_start_struct = cell2struct(num2cell(b.X_starts),b.var_names(1:end-1)',1);
% struct('D_f',b.D_f_start,'D_s_ratio',b.D_s_ratio_start,...
%         'h_f_ratio',b.h_f_ratio_start,'T_s_ratio',b.T_s_ratio_start,...
%         'F_max',b.F_max_start,'D_int',b.B_p_start,'w_n',b.w_n_start);

b.constraint_names = {'float_too_heavy','float_too_light','spar_too_heavy','spar_too_light',...
                      'stability','FOS_float_yield','FOS_col_yield','FOS_plate','FOS_col_buckling',...
                      'pos_power','spar_damping','spar_height_up','spar_height_down','LCOE_max',...
                      'irrelevant_max_force','water_deep_enough'};

% modify nominal control inputs to so power and force matches actual
[F_max_nom, B_p_nom, w_n_nom] = deal(1e6);%find_nominal_inputs(b, false);
b.F_max_nom = F_max_nom;
b.B_p_nom = B_p_nom;
b.w_n_nom = w_n_nom;

b.X_noms = [b.D_f_nom b.D_s_ratio_nom b.h_f_ratio_nom b.T_s_ratio_nom b.F_max_nom b.B_p_nom b.w_n_nom]';

end