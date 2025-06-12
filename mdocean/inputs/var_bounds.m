function b = var_bounds(mode)
% mode = 'wecsim': use parameters corresponding to RM3.out in WEC-Sim
% mode = anything else or not provided: use parmaters corresponding to RM3 report (default)

if nargin<1
    mode = '';
end

% design variables
b.var_names = {'D_s','D_f','T_f_2','h_s','h_fs_clear','F_max','P_max',...
                't_fb','t_sr','t_d','h_stiff_f','h1_stiff_d','M'};
b.var_names_pretty = {'D_s','D_f','T_{f,2}','h_s','h_{fs,clear}','F_{max}','P_{max}',...
                't_{fb}','t_{sr}','t_d','h_{stiff,f}','h_{1,stiff,d}','M'};
b.var_descs = {'Spar Diameter','Float Diameter','Float Draft','Spar Height',...
    'Float-Spar Clearance Height','Maximum Force','Maximum Power',...
    'Float Bottom Thickness','Spar Radial Thickness','Damping Plate Thickness',...
    'Float Stiffener Height','Damping Plate Stiffener Height','Material Index'};

% diameter of spar (m)	
b.D_s_min = 0;
b.D_s_max = 30;
b.D_s_nom = 6;
b.D_s_wecsim = 6;
b.D_s_start = 6;

% outer diameter of float (m)	
b.D_f_min = 1;
b.D_f_max = 80;
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

% vertical clearance between float tube support and spar (m)
b.h_fs_clear_min = .01;
b.h_fs_clear_max = 10;
b.h_fs_clear_wecsim = 4;
b.h_fs_clear_nom = 4; % p169 11/26/24
b.h_fs_clear_start = 4;

% maximum powertrain force (MN)
b.F_max_min = 0.01;
b.F_max_max = 100;
b.F_max_wecsim = 100;
b.F_max_nom = 1;
b.F_max_start = 5;

% maximum generator power (100 kW)
b.P_max_min = .01;
b.P_max_max = 30;
b.P_max_nom = 2.86;
b.P_max_wecsim = 2.86;
b.P_max_start = 2.86;

in2mm = 25.4;
% material thickness of float bottom (mm)
b.t_fb_min = 0.05 * in2mm;
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

% material thickness of damping plate (mm)
b.t_d_min = 0.05 * in2mm;
b.t_d_max = 2.0 * in2mm;
b.t_d_nom = 1.0 * in2mm;
b.t_d_wecsim = 1.0 * in2mm;
b.t_d_start = 1.0 * in2mm;

in2m = in2mm / 1000;
% float stiffener height (dm)
b.h_stiff_f_min = 0;
b.h_stiff_f_max = 3*10;
b.h_stiff_f_nom    = 16 * in2m * 10;
b.h_stiff_f_wecsim = 16 * in2m * 10;
b.h_stiff_f_start  = 16 * in2m * 10;

% damping plate stiffener height (m)
b.h1_stiff_d_min = 0;
b.h1_stiff_d_max = 2;
b.h1_stiff_d_nom    = 22 * in2m;
b.h1_stiff_d_wecsim = 22 * in2m;
b.h1_stiff_d_start  = 22 * in2m;

% material index (-)
b.M_min = 1;
b.M_max = 3;
b.M_wecsim = 1;
b.M_nom = 1;
b.M_start = 1;

                 % D_s    D_f   T_f_2   h_s    h_fs_clear  F_max  P_max  thicknesses  stiffener-heights ]
b.mins_flexible = [false  true  true    true   false       true   true   true(1,3)    false(1,2)]';
b.maxs_flexible = [true   true  false   false  true        true   true   true(1,3)    true(1,2) ]';
% if a bound is marked flexible and the bound is active after optimization, 
% a warning in gradient_optim will remind you to adjust the bound.

% create vectors of mins, maxs, noms, and starts
n_dv = length(b.var_names)-1;
[X_mins,X_maxs,X_noms,X_noms_wecsim,X_starts] = deal(zeros(n_dv,1));
for i=1:n_dv
    dv_name = b.var_names{i};
    X_mins(i) = b.([dv_name '_min']);
    X_maxs(i) = b.([dv_name '_max']);
    X_noms(i) = b.([dv_name '_nom']);
    X_noms_wecsim(i) = b.([dv_name '_wecsim']);
    X_starts(i) = b.([dv_name '_start']);
end
b.X_mins = X_mins;
b.X_maxs = X_maxs;
if strcmpi(mode,'wecsim')
    b.X_noms = X_noms_wecsim;
else
    b.X_noms = X_noms;
end
b.X_starts = X_starts;

b.X_start_struct = cell2struct(num2cell(b.X_starts),b.var_names(1:end-1)',1);

% constraints
b.constraint_names = {'float_too_heavy','float_too_light','spar_too_heavy','spar_too_light',...
                      'stability','FOS_float_max','FOS_float_fatigue',...
                      'FOS_col_max','FOS_col_fatigue','FOS_plate_max','FOS_plate_fatigue',...
                      'FOS_col_local_max','FOS_col_local_fatigue',...
                      'pos_power','LCOE_max','irrelevant_max_force',...
                      'force_limit','power_limit',...
                      'spar_height_up','spar_height_down','float_spar_hit',...
                      'linear_theory'};
i1 = length(b.constraint_names);
JPD_size = 14*15;
storm_size = 7;
for i = (i1+1):(i1+JPD_size+storm_size)
    b.constraint_names{i} = strcat('prevent_slamming',num2str(i-i1));
end
b.constraint_names_pretty = remove_underscores(b.constraint_names);

b.lin_constraint_names = {'spar_natural_freq','float_spar_diam','float_spar_draft',...
                          'float_spar_tops','float_seafloor','spar_seafloor',...
                          'damping_plate_thickness','float_stiffener'};
b.lin_constraint_names_pretty = remove_underscores(b.lin_constraint_names);

% objectives
b.obj_names = {'LCOE','capex_design'};
b.obj_names_pretty = {'LCOE','C_{design}'};

% indices
[~,idxs_sort] = sort(b.var_names(1:end-1)); % alphabetical design variable indices
idxs_recover = zeros(size(idxs_sort));
idxs_recover(idxs_sort) = 1:length(idxs_sort); % indices to recover unsorted variabes from sorted ones
b.idxs_sort    = idxs_sort;
b.idxs_recover = idxs_recover;

% uuid
b.filename_uuid = ''; % string to append to generated filenames to prevent parallel overlap

% calibrations of nominal values
b.F_max_nom = find_nominal_inputs(b, parameters(mode));
b.X_noms(strcmp(b.var_names,'F_max')) = b.F_max_nom;

% balanced design
b.power_balanced = 140e3;

end