function b = var_bounds()

b.var_names = {'D_s','D_f','T_f_2','h_s','F_max','B_p','w_n',...
                't_ft','t_fr','t_fc','t_fb','t_sr','t_dt','P_max','M'};
b.var_names_pretty = {'D_s','D_f','T_{f,2}','h_s','F_{max}','B_p','\omega_n',...
                't_{ft}','t_{fr}','t_{fc}','t_{fb}','t_{sr}','t_{dt}','P_{max}','M'};

% inner diameter of float (m)	
b.D_s_min = 0;
b.D_s_max = 30;
b.D_s_nom = 6;
b.D_s_start = 6;

% outer diameter of float (m)	
b.D_f_min = 1;
b.D_f_max = 40;
b.D_f_nom = 20;
b.D_f_start = 20;

% draft of float (m)
b.T_f_2_min = .5;
b.T_f_2_max = 100;
b.T_f_2_nom = 3;
b.T_f_2_start = 3;

% height of spar (m)
b.h_s_min = 5;
b.h_s_max = 100;
b.h_s_nom = 37.9;%44;
b.h_s_start = 37.9;%44;

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

% material thickness of float top (m)
b.t_ft_min = 0.3;
b.t_ft_max = 0.5;
b.t_ft_nom = 0.5;
b.t_ft_start = 0.4;

% material thickness of float bottom (m)
b.t_fr_min = 0.1; 
b.t_fr_max = 0.6;
b.t_fr_nom = 0.44;
b.t_fr_start = 0.2;

% materal thickness of float circumferential gussets (m)
b.t_fc_min = 0.1; 
b.t_fc_max = 0.6; 
b.t_fc_nom = 0.44; 
b.t_fc_start = 0.2; 

% material thickness of float bottom (m)
b.t_fb_min = 0.1;
b.t_fb_max = 0.9;
b.t_fb_nom = 0.56;
b.t_fb_start = 0.2;

% material thickness of spar radial (m)
b.t_sr_min = 0.2;
b.t_sr_max = 2.0;
b.t_sr_nom = 1.00;
b.t_sr_start = 0.3;

% material thickness of damping plate (m)
b.t_dt_min = 0.2;
b.t_dt_max = 2.0;
b.t_dt_nom = 1.00;
b.t_dt_start = 0.3;

% maximum generator power (W)
b.P_max_min = 3;
b.P_max_max = 40;
b.P_max_nom = 20;
b.P_max_start = 5;

% material index (-)
b.M_min = 1;
b.M_max = 3;
b.M_nom = 1;
b.M_start = 1;

                 % D_s    D_f   T_f_2  h_s    F_max  B_p   w_n  thicknesses P_max]
b.mins_flexible = [false  true  true   true   true   true  true true(6,1)   true]';
b.maxs_flexible = [true   true  false  false  true   true  true true(6,1)   true]';
% if a bound is marked flexible and the bound is active after optimization, 
% a warning in gradient_optim will remind you to adjust the bound.

b.X_mins = [b.D_s_min b.D_f_min b.T_f_2_min b.h_s_min b.F_max_min b.B_p_min b.w_n_min ...
    b.t_ft_min b.t_fr_min b.t_fc_min b.t_fb_min b.t_sr_min b.t_dt_min b.P_max_min]';
b.X_maxs = [b.D_s_max b.D_f_max b.T_f_2_max b.h_s_max b.F_max_max b.B_p_max b.w_n_max ...
    b.t_ft_max b.t_fr_max b.t_fc_max b.t_fb_max b.t_sr_max b.t_dt_max b.P_max_max]';
b.X_noms = [b.D_s_nom b.D_f_nom b.T_f_2_nom b.h_s_nom b.F_max_nom b.B_p_nom b.w_n_nom ...
    b.t_ft_nom b.t_fr_nom b.t_fc_nom b.t_fb_nom b.t_sr_nom b.t_dt_nom b.P_max_nom]';
b.X_starts = [b.D_s_start b.D_f_start b.T_f_2_start b.h_s_start b.F_max_start b.B_p_start b.w_n_start ...
    b.t_ft_start b.t_fr_start b.t_fc_start b.t_fb_start b.t_sr_start b.t_dt_start b.P_max_start]';

b.X_start_struct = cell2struct(num2cell(b.X_starts),b.var_names(1:end-1)',1);

b.constraint_names = {'float_too_heavy','float_too_light','spar_too_heavy','spar_too_light',...
                      'stability','FOS_float_yield','FOS_col_yield','FOS_plate','FOS_col_buckling',...
                      'pos_power','LCOE_max','irrelevant_max_force',...
                      'spar_height_up','spar_height_down','linear_theory'};
i1 = length(b.constraint_names);
for i = (i1+1):(i1+14*15)
    b.constraint_names{i} = strcat('prevent_slamming',num2str(i-i1));
end

b.lin_constraint_names = {'spar_natural_freq','float_spar_diam','float_spar_draft',
                          'float_spar_tops','float_seafloor','spar_seafloor'};

[~,idxs_sort] = sort(b.var_names(1:end-1)); % alphabetical design variable indices
idxs_recover = zeros(size(idxs_sort));
idxs_recover(idxs_sort) = 1:length(idxs_sort); % indices to recover unsorted variabes from sorted ones
b.idxs_sort    = idxs_sort;
b.idxs_recover = idxs_recover;

b.filename_uuid = ''; % string to append to generated filenames to prevent parallel overlap

end