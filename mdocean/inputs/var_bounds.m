function b = var_bounds(mode)
% mode = 'wecsim': use parameters corresponding to RM3.out in WEC-Sim
% mode = anything else or not provided: use parmaters corresponding to RM3 report (default)

if nargin<1
    mode = '';
end

b.var_names = {'D_s','D_f','T_f_2','h_s','F_max','B_p','w_n','M'};
b.var_names_pretty = {'D_s','D_f','T_{f,2}','h_s','F_{max}','B_p','\omega_n','M'};

% inner diameter of float (m)	
b.D_s_min = 0;
b.D_s_max = 30;
b.D_s_nom = 6.5;
b.D_s_wecsim = 6;
b.D_s_start = 6;

% outer diameter of float (m)	
b.D_f_min = 1;
b.D_f_max = 40;
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
b.F_max_wecsim = 1e6;
b.F_max_nom = 100;
b.F_max_start = 5;

% powertrain damping (MN / (m/s))
b.B_p_min = .1;
b.B_p_max = 50;
b.B_p_wecsim = 10;
b.B_p_nom = 10;
b.B_p_start = 0.5;

% natural frequency (rad/s)
b.w_n_min = .01;%2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy
b.w_n_max = 40;%2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
b.w_n_wecsim = 0.8;
b.w_n_nom = 0.8;
b.w_n_start = 0.8;

% material index (-)
b.M_min = 1;
b.M_max = 3;
b.M_wecsim = 1;
b.M_nom = 1;
b.M_start = 1;

                 % D_s    D_f   T_f_2  h_s    F_max  B_p   w_n]
b.mins_flexible = [false  true  true   true   true   true  true]';
b.maxs_flexible = [true   true  false  false  true   true  true]';
% if a bound is marked flexible and the bound is active after optimization, 
% a warning in gradient_optim will remind you to adjust the bound.

b.X_mins = [b.D_s_min b.D_f_min b.T_f_2_min b.h_s_min b.F_max_min b.B_p_min b.w_n_min]';
b.X_maxs = [b.D_s_max b.D_f_max b.T_f_2_max b.h_s_max b.F_max_max b.B_p_max b.w_n_max]';
if strcmpi(mode,'wecsim')
    b.X_noms = [b.D_s_wecsim b.D_f_wecsim b.T_f_2_wecsim b.h_s_wecsim b.F_max_wecsim b.B_p_wecsim b.w_n_wecsim]';
else
    b.X_noms = [b.D_s_nom b.D_f_nom b.T_f_2_nom b.h_s_nom b.F_max_nom b.B_p_nom b.w_n_nom]';
end
b.X_starts = [b.D_s_start b.D_f_start b.T_f_2_start b.h_s_start b.F_max_start b.B_p_start b.w_n_start]';

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