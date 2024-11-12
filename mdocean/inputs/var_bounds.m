function b = var_bounds()

b.var_names = {'D_s','D_s_over_D_f','T_f_over_T_s','T_s_over_h_s','F_max','B_p','w_n','M'};
b.var_names_pretty = {'D_s','D_s/D_f','T_f/T_s','T_s/h_s','F_{max}','B_p','\omega_n','M'};

% outer diameter of float (m)	
b.D_s_min = 6;    % a consequence of the spar natural frequency constraint
b.D_s_max = 30;
b.D_s_nom = 6;
b.D_s_start = 6;

% inner diameter to outer diameter of float ratio (-)	
b.D_s_over_D_f_min = 0.01;
b.D_s_over_D_f_max = 0.99;
b.D_s_over_D_f_nom = 6/20;
b.D_s_over_D_f_start = 6/20;

% height of float to outer diameter of float ratio (-)
b.T_f_over_T_s_min = .01;
b.T_f_over_T_s_max = 0.99;
b.T_f_over_T_s_nom = 3/29;
b.T_f_over_T_s_start = 3/29;

% draft of spar to height of spar ratio (-)
b.T_s_over_h_s_min = 0.01;
b.T_s_over_h_s_max = 0.99;
b.T_s_over_h_s_nom = 29/37.9;%35/44;
b.T_s_over_h_s_start = 29/37.9;%35/44;

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

b.X_mins = [b.D_s_min b.D_s_over_D_f_min b.T_f_over_T_s_min b.T_s_over_h_s_min b.F_max_min b.B_p_min b.w_n_min]';
b.X_maxs = [b.D_s_max b.D_s_over_D_f_max b.T_f_over_T_s_max b.T_s_over_h_s_max b.F_max_max b.B_p_max b.w_n_max]';
b.X_noms = [b.D_s_nom b.D_s_over_D_f_nom b.T_f_over_T_s_nom b.T_s_over_h_s_nom b.F_max_nom b.B_p_nom b.w_n_nom]';
b.X_starts = [b.D_s_start b.D_s_over_D_f_start b.T_f_over_T_s_start b.T_s_over_h_s_start b.F_max_start b.B_p_start b.w_n_start]';

b.X_start_struct = cell2struct(num2cell(b.X_starts),b.var_names(1:end-1)',1);

b.constraint_names = {'float_too_heavy','float_too_light','spar_too_heavy','spar_too_light',...
                      'stability','FOS_float_yield','FOS_col_yield','FOS_plate','FOS_col_buckling',...
                      'pos_power','spar_damping','LCOE_max','irrelevant_max_force','water_deep_enough',...
                      'spar_height_up','spar_height_down','linear_theory'};
for i = 18:(17+14*15)
    b.constraint_names{i} = strcat('prevent_slamming',num2str(i-17));
end

[~,idxs_sort] = sort(b.var_names(1:end-1)); % alphabetical design variable indices
idxs_recover = zeros(size(idxs_sort));
idxs_recover(idxs_sort) = 1:length(idxs_sort); % indices to recover unsorted variabes from sorted ones
b.idxs_sort    = idxs_sort;
b.idxs_recover = idxs_recover;

b.filename_uuid = ''; % string to append to generated filenames to prevent parallel overlap

end