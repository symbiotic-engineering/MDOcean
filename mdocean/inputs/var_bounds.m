function b = var_bounds(mode)
% mode = 'wecsim': use parameters corresponding to RM3.out in WEC-Sim
% mode = anything else or not provided: use parmaters corresponding to RM3 report (default)

if nargin<1
    mode = '';
end

b.var_names = {'D_s','D_f','T_f_2','h_s','F_max','B_p','h_fs_clear',...
                't_ft','t_fr','t_fc','t_fb','t_sr','t_dt','P_max','M'};
b.var_names_pretty = {'D_s','D_f','T_{f,2}','h_s','F_{max}','B_p','h_{fs,clear}',...
                't_{ft}','t_{fr}','t_{fc}','t_{fb}','t_{sr}','t_{dt}','P_{max}','M'};

% inner diameter of float (m)	
b.D_s_min = 0;
b.D_s_max = 30;
b.D_s_nom = 6;
b.D_s_wecsim = 6;
b.D_s_start = 6;

% outer diameter of float (m)	
b.D_f_min = 1;
b.D_f_max = 50;
b.D_f_wecsim = 20;
b.D_f_nom = 20;
b.D_f_start = 20;

% draft of float (m)
b.T_f_2_min = .5;
b.T_f_2_max = 100;
b.T_f_2_nom = 3.2;
b.T_f_2_wecsim = 3;
b.T_f_2_start = 3;

% height of spar (m)
b.h_s_min = 5;
b.h_s_max = 100;
b.h_s_nom = 44;
b.h_s_wecsim = 37.9;
b.h_s_start = 37.9;

% maximum powertrain force (MN)
b.F_max_min = 0.01;
b.F_max_max = 100;
b.F_max_wecsim = 100;
b.F_max_nom = 100;
b.F_max_start = 5;

% powertrain damping (MN / (m/s))
b.B_p_min = .1;
b.B_p_max = 50;
b.B_p_wecsim = 10;
b.B_p_nom = 10;
b.B_p_start = 0.5;

% vertical clearance between float tube support and spar (m)
b.h_fs_clear_min = .01;
b.h_fs_clear_max = 10;
b.h_fs_clear_wecsim = 4;
b.h_fs_clear_nom = 4; % p169 11/26/24
b.h_fs_clear_start = 4;

in2mm = 25.4;
% material thickness of float top (mm)
b.t_ft_min = 0.1 * in2mm;
b.t_ft_max = 1.0 * in2mm;
b.t_ft_nom = 0.5 * in2mm;
b.t_ft_wecsim = 0.5 * in2mm;
b.t_ft_start = 0.5 * in2mm;

% material thickness of float bottom (mm)
b.t_fr_min = 0.1 * in2mm; 
b.t_fr_max = 1.0 * in2mm;
b.t_fr_nom = 0.44 * in2mm;
b.t_fr_wecsim = 0.44 * in2mm;
b.t_fr_start = 0.44 * in2mm;

% materal thickness of float circumferential gussets (mm)
b.t_fc_min = 0.1 * in2mm; 
b.t_fc_max = 1.0 * in2mm; 
b.t_fc_nom = 0.44 * in2mm; 
b.t_fc_wecsim = 0.44 * in2mm; 
b.t_fc_start = 0.44 * in2mm; 

% material thickness of float bottom (mm)
b.t_fb_min = 0.1 * in2mm;
b.t_fb_max = 1.0 * in2mm;
b.t_fb_nom = 0.56 * in2mm;
b.t_fb_wecsim = 0.56 * in2mm;
b.t_fb_start = 0.56 * in2mm;

% material thickness of spar radial (mm)
b.t_sr_min = 0.2 * in2mm;
b.t_sr_max = 2.0 * in2mm;
b.t_sr_nom = 1.0 * in2mm;
b.t_sr_wecsim = 1.0 * in2mm;
b.t_sr_start = 1.0 * in2mm;

% material thickness of damping plate support tube radial walls (mm)
b.t_dt_min = 0.2 * in2mm;
b.t_dt_max = 2.0 * in2mm;
b.t_dt_nom = 1.0 * in2mm;
b.t_dt_wecsim = 1.0 * in2mm;
b.t_dt_start = 1.0 * in2mm;

% maximum generator power (kW)
b.P_max_min = 50;
b.P_max_max = 1000;
b.P_max_nom = 286;
b.P_max_wecsim = 286;
b.P_max_start = 286;

% material index (-)
b.M_min = 1;
b.M_max = 3;
b.M_wecsim = 1;
b.M_nom = 1;
b.M_start = 1;

                 % D_s    D_f   T_f_2  h_s    F_max  B_p   h_fs_clear  thicknesses P_max]
b.mins_flexible = [false  true  true   true   true   true  true        true(1,6)   true]';
b.maxs_flexible = [true   true  false  false  true   true  true        true(1,6)   true]';
% if a bound is marked flexible and the bound is active after optimization, 
% a warning in gradient_optim will remind you to adjust the bound.

b.X_mins = [b.D_s_min b.D_f_min b.T_f_2_min b.h_s_min b.F_max_min b.B_p_min b.h_fs_clear_min ...
    b.t_ft_min b.t_fr_min b.t_fc_min b.t_fb_min b.t_sr_min b.t_dt_min b.P_max_min]';
b.X_maxs = [b.D_s_max b.D_f_max b.T_f_2_max b.h_s_max b.F_max_max b.B_p_max b.h_fs_clear_max ...
    b.t_ft_max b.t_fr_max b.t_fc_max b.t_fb_max b.t_sr_max b.t_dt_max b.P_max_max]';

if strcmpi(mode,'wecsim')
    b.X_noms = [b.D_s_wecsim b.D_f_wecsim b.T_f_2_wecsim b.h_s_wecsim b.F_max_wecsim ...
        b.B_p_wecsim b.h_fs_clear_wecsim b.t_ft_wecsim b.t_fr_wecsim b.t_fc_wecsim ...
        b.t_fb_wecsim b.t_sr_wecsim b.t_dt_wecsim b.P_max_wecsim]';
else
    b.X_noms = [b.D_s_nom b.D_f_nom b.T_f_2_nom b.h_s_nom b.F_max_nom b.B_p_nom b.h_fs_clear_nom ...
        b.t_ft_nom b.t_fr_nom b.t_fc_nom b.t_fb_nom b.t_sr_nom b.t_dt_nom b.P_max_nom]';
end

b.X_starts = [b.D_s_start b.D_f_start b.T_f_2_start b.h_s_start b.F_max_start b.B_p_start b.h_fs_clear_start ...
    b.t_ft_start b.t_fr_start b.t_fc_start b.t_fb_start b.t_sr_start b.t_dt_start b.P_max_start]';

b.X_start_struct = cell2struct(num2cell(b.X_starts),b.var_names(1:end-1)',1);

b.constraint_names = {'float_too_heavy','float_too_light','spar_too_heavy','spar_too_light',...
                      'stability','FOS_float_max','FOS_float_fatigue',...
                      'float_top_ratio','float_radial_ratio','float_circ_ratio',... % placeholders for now
                      'FOS_col_max','FOS_col_fatigue','FOS_plate_max','FOS_plate_fatigue',...
                      'FOS_col_local_max','FOS_col_local_fatigue',...
                      'pos_power','LCOE_max','irrelevant_max_force',...
                      'spar_height_up','spar_height_down','float_spar_hit',...
                      'linear_theory'};
i1 = length(b.constraint_names);
for i = (i1+1):(i1+14*15)
    b.constraint_names{i} = strcat('prevent_slamming',num2str(i-i1));
end
b.constraint_names_pretty = remove_underscores(b.constraint_names);

b.lin_constraint_names = {'spar_natural_freq','float_spar_diam','float_spar_draft',...
                          'float_spar_tops','float_seafloor','spar_seafloor'};
b.lin_constraint_names_pretty = remove_underscores(b.lin_constraint_names);

[~,idxs_sort] = sort(b.var_names(1:end-1)); % alphabetical design variable indices
idxs_recover = zeros(size(idxs_sort));
idxs_recover(idxs_sort) = 1:length(idxs_sort); % indices to recover unsorted variabes from sorted ones
b.idxs_sort    = idxs_sort;
b.idxs_recover = idxs_recover;

b.filename_uuid = ''; % string to append to generated filenames to prevent parallel overlap

b.F_max_nom = find_nominal_inputs(b, mode, false);
b.X_noms(5) = b.F_max_nom;

b.obj_names = {'LCOE','capex_design'};
b.obj_names_pretty = {'LCOE','C_{design}'};

end