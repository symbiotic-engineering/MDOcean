function [J, P_matrix_elec, g, val] = simulation(X, p)	

X = max(X,1e-3); % sometimes optimizer will choose inputs that violate bounds, this is to prevent errors from negative numbers

%% Assemble inputs
in = p;

in.D_s   = X(1);        % inner diameter of float (m)
in.D_f   = X(2);        % outer diameter of float (m)
in.T_f_2 = X(3);        % draft of float (m)
in.h_s   = X(4);        % total height of spar (m)
in.h_fs_clear = X(5);   % vertical clearance between float tubes and spar when at rest (m)
in.F_max = X(6) * 1e6;  % max powertrain force (N)
in.P_max = X(7) * 1e5;  % maximum power (W)
in.t_f_b = X(8) * 1e-3; % float bottom thickness (m)
in.t_s_r = X(9) * 1e-3; % vertical column thickness (m)
in.t_d   = X(10)* 1e-3; % damping plate thickness (m)
in.h_stiff_f  = X(11)/10;  % float stiffener height (m)
in.h1_stiff_d = X(12);  % damping plate stiffener larger height (m)
in.M     = X(13);       % material index (-)

% float thickness ratios
in.t_f_r  = in.t_f_b * p.t_f_r_over_t_f_b;  % float radial wall thickness (m)
in.t_f_c  = in.t_f_b * p.t_f_c_over_t_f_b;  % float circumferential gusset thickness (m)
in.t_f_t  = in.t_f_b * p.t_f_t_over_t_f_b;  % float top thickness (m)

% stiffener thickness ratios
in.w_stiff_f = p.w_over_h_stiff_f  * in.h_stiff_f;
in.h_stiff_d = p.h_over_h1_stiff_d * in.h1_stiff_d;
in.w_stiff_d = p.w_over_h1_stiff_d * in.h1_stiff_d;

% Geometric similarity to maintain constant damping ratio
% D_s sets D_d, T_s, h_d
in.D_d = p.D_d_over_D_s * in.D_s;
in.T_s = p.T_s_over_D_s * in.D_s;
in.h_d = p.h_d_over_D_s * in.D_s;

% Geometric similarity to float submergence parameter
in.h_f = in.T_f_2 / p.T_f_2_over_h_f;

% Geometric similarity to float angle paramter
in.T_f_1 = p.T_f_1_over_T_f_2 * in.T_f_2;
in.D_f_b = p.D_f_b_over_D_f * in.D_f;

% Space between float and spar
D_f_in = p.D_f_in_over_D_s * in.D_s;
in.D_f_in = min(D_f_in, in.D_f - 1e-6); % ensure positive area to avoid sim breaking for infeasible inputs

%% Run modules
[V_d, m_m, m_f_tot, m_s_tot, ...
 A_c, A_lat_sub, ...
 I, T, V_f_pct, V_s_pct, GM, ...
 A_dt, L_dt, A_hull] = geometry(in.D_s, in.D_f, in.D_f_in, in.D_f_b, ...
                                in.T_f_1, in.T_f_2, in.h_f, in.h_s, ...
                                in.h_fs_clear, in.D_f_tu, in.t_f_t, ...
                                in.t_f_r, in.t_f_c, in.t_f_b, in.t_f_tu, ...
                                in.t_s_r, in.t_d_tu, in.D_d, in.D_d_tu, ...
                                in.theta_d_tu, in.T_s, in.h_d, in.t_d,...
                                in.h_stiff_f, in.w_stiff_f, in.num_sections_f, ...
                                in.h_stiff_d, in.w_stiff_d, in.num_stiff_d, ...
                                in.M, in.rho_m, in.rho_w, in.m_scale);

m_f_tot = max(m_f_tot,1e-3); % zero out negative mass produced by infeasible inputs

[F_heave_storm, F_surge_storm, ...
 F_heave_op, F_surge_op, F_ptrain_max, ...
 P_var, P_avg_elec, P_matrix_elec, ...
 X_constraints,...
 B_p,K_p,mag_U,X_u,X_f,X_s,...
 P_matrix_mech,P_sat_ratio,...
 A_f_over_rho, A_s_over_rho,A_c_over_rho,...
 B_f_over_rho_w, B_s_over_rho_w,B_c_over_rho_w,...
 gamma_f_over_rho_g,gamma_s_over_rho_g,...
 gamma_phase_f,gamma_phase_s,w,...
 phase_X_f,phase_X_u]           = dynamics(in, m_f_tot, m_s_tot);

[FOS_float,FOS_spar,FOS_damping_plate,...
    FOS_spar_local] = structures(...
          	F_heave_storm, F_surge_storm, F_heave_op, F_surge_op, ... % forces
            in.h_s, in.T_s, in.D_s, in.D_f, in.D_f_in, in.num_sections_f, ... % bulk dimensions
            in.D_f_tu, in.D_d, L_dt, in.theta_d_tu, in.D_d_tu,... % more bulk dimensions
            in.t_s_r, I, A_c, A_lat_sub, in.t_f_b, in.t_f_t, in.t_d, in.t_d_tu, in.h_d, A_dt,  ... % structural dimensions
            in.h_stiff_f, in.w_stiff_f, in.h_stiff_d, in.w_stiff_d,... % stiffener thicknesses
            in.M, in.rho_w, in.g, in.sigma_y, in.sigma_e, in.E, in.nu, ... % constants
            in.num_terms_plate, in.radial_mesh_plate, in.num_stiff_d);

[LCOE,capex_design,capex] = econ(m_m, in.M, in.cost_perkg_mult, in.N_WEC, P_avg_elec, ...
                            in.FCR, in.cost_perN_mult, in.cost_perW_mult, ...
                            in.F_max, in.P_max, in.eff_array);

model_grid_value = false; % override for merging
if model_grid_value
    [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = grid_CEM(B_p, X_u, phase_X_u, ...
                                                    gamma_phase_f, capex/in.P_max*1000, in.location);
    
    net_eco_value = enviro(m_m, in.distance_from_shore, A_hull, ...
                            CEM_CO2, CEM_wec_capacity, in.P_max, ...
                            in.eco_cost_steel, in.eco_cost_fiberglass, ...
                            in.eco_cost_distance, in.eco_cost_carbon);
else
    [CEM_CO2, CEM_wec_capacity, CEM_grid_cost,net_eco_value] = deal(0);
end


J_capex_design = capex_design / 1e6; % convert $ to $M

J = [LCOE, J_capex_design, CEM_grid_cost, net_eco_value];

%% Assemble constraints g(x) >= 0
num_g = 23+numel(p.JPD)+length(p.T_struct);
g = zeros(1,num_g);
g(1) = V_f_pct;                         % prevent float too heavy
g(2) = 1 - V_f_pct;                     % prevent float too light
g(3) = V_s_pct;                         % prevent spar too heavy
g(4) = 1 - V_s_pct;                     % prevent spar too light
g(5) = GM;                              % pitch stability of float-spar system
g(6) = FOS_float(1) / p.FOS_min - 1;    % float survives max force
g(7) = FOS_float(2) / p.FOS_min - 1;    % float survives fatigue
g(8) = FOS_spar(1) / p.FOS_min - 1;     % spar survives max force
g(9) = FOS_spar(2) / p.FOS_min - 1;     % spar survives fatigue
g(10) = FOS_damping_plate(1) * in.FOS_mult_d / p.FOS_min - 1; % damping plate survives max force
g(11) = FOS_damping_plate(2) * in.FOS_mult_d / p.FOS_min - 1; % damping plate survives fatigue
g(12) = FOS_spar_local(1) / p.FOS_min - 1;    % spar survives max force in local buckling
g(13) = FOS_spar_local(2) / p.FOS_min - 1;    % spar survives fatigue in local buckling
if ~isfield(p,'avg_power_min') || p.avg_power_min == 0
    g(14) = P_avg_elec/1e6;                     % positive power
else
    g(14) = P_avg_elec/p.avg_power_min - 1; % prevent less avg power than threshold
end
%1 + min(Kp_over_Ks,[],'all');   % spar heave stability (positive effective stiffness)
g(15) = p.LCOE_max/LCOE - 1;            % prevent more expensive than threshold
g(16) = F_ptrain_max/in.F_max - 1;      % prevent irrelevant max force -
                                        % this constraint should always be
                                        % active unless F_max==Inf
                                        % and is only required when p.cost_perN = 0.
g(17) = in.F_max/F_ptrain_max - 1;      % prevent more PTO force than is available
                                        % this constraint should always be
                                        % active unless F_max==Inf and is
                                        % only required when p.use_force_sat = false.
g(18) = in.P_max/max(P_matrix_elec,[],'all') - 1; % prevent more PTO power than is avaliable
                                        % this constraint should always be active
                                        % and is only required when p.use_power_sat = false.
g(19) = X_constraints(1);               % prevent float rising above top of spar
g(20) = X_constraints(2);               % prevent float going below bottom of spar
g(21) = X_constraints(3);               % prevent float support tube (PTO attachment) from hitting spar
g(22) = X_constraints(4);               % float amplitude obeys linear theory
g(23) = X_constraints(5);               % spar amplitude obeys linear theory
g(24:end) = X_constraints(6:end);       % prevent rising out of water/slamming

no_inf = all(~isinf([g LCOE P_var]));
no_nan = all(~isnan([g LCOE P_var]));
no_imag = all(isreal([g LCOE P_var]));
criteria = no_inf && no_nan && no_imag;
if ~criteria
    disp(['g: ' num2str(g)])
    disp(['LCOE: ' num2str(LCOE) ', P_var: ' num2str(P_var)])
    warning('Inf, NaN, or imaginary constraint or objective detected')
end

if nargout > 3 % if returning extra struct output for validation
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,~,...
     mass,CB_f,CG_f] = geometry(in.D_s, in.D_f, in.D_f_in, in.D_f_b, ...
                                 in.T_f_1, in.T_f_2, in.h_f, in.h_s, ...
                                 in.h_fs_clear, in.D_f_tu, in.t_f_t, ...
                                 in.t_f_r, in.t_f_c, in.t_f_b, in.t_f_tu, ...
                                 in.t_s_r, in.t_d_tu, in.D_d, in.D_d_tu, ...
                                 in.theta_d_tu, in.T_s, in.h_d, in.t_d,...
                                 in.h_stiff_f, in.w_stiff_f, in.num_sections_f, ...
                                 in.h_stiff_d, in.w_stiff_d, in.num_stiff_d, ...
                                 in.M, in.rho_m, in.rho_w, in.m_scale);

    [~,~,capex,opex,pto, devicestructure] = econ(m_m, in.M, in.cost_perkg_mult, in.N_WEC, P_avg_elec, in.FCR, ...
                        in.cost_perN_mult, in.cost_perW_mult, in.F_max, in.P_max, in.eff_array);

    [~, ~, ~, ~, ~, ~, ~, ~, ~, ...
     B_p,K_p,mag_U,X_u,X_f,X_s,...
     P_matrix_mech,P_sat_ratio,...
     A_f_over_rho, A_s_over_rho,A_c_over_rho,...
     B_f_over_rho_w, B_s_over_rho_w,B_c_over_rho_w,...
     gamma_f_over_rho_g,gamma_s_over_rho_g,...
     gamma_phase_f,gamma_phase_s,w,...
     phase_X_f,phase_X_u] = dynamics(in, m_f_tot, m_s_tot);

    val.mass_f  = mass(1);
    val.mass_vc = mass(2);
    val.mass_rp = mass(3);
    val.mass_tot = sum(mass);
    val.capex = capex;
    val.capex_design = capex_design;
    val.J_capex_design = J_capex_design;
    val.capex_struct = devicestructure;
    val.capex_PTO = pto;
    val.opex = opex;
    val.LCOE = LCOE;
    val.power_avg = P_avg_elec;
    val.power_max = max(P_matrix_elec,[],'all');
    val.force_heave = F_heave_storm;
    val.force_heave_op = F_heave_op;
    val.force_ptrain = F_ptrain_max;
    val.FOS_spar = FOS_spar(1);
	val.c_v = P_var;
    val.B_p = B_p;
    val.K_p = K_p;
    val.mag_U = mag_U;
    val.X_u = X_u;
    val.X_f = X_f;
    val.X_s = X_s;
    val.P_mech = P_matrix_mech;
    val.P_sat_ratio = P_sat_ratio;
    val.CB_f = CB_f;
    val.CG_f = CG_f;
    val.vol_f = V_d(1);
    val.vol_s = V_d(2) + V_d(3);
    val.A_f_over_rho = A_f_over_rho;
    val.A_s_over_rho = A_s_over_rho;
    val.A_c_over_rho = A_c_over_rho;
    val.B_f_over_rho_w = B_f_over_rho_w;
    val.B_s_over_rho_w = B_s_over_rho_w;
    val.B_c_over_rho_w = B_c_over_rho_w;
    val.gamma_f_over_rho_g = gamma_f_over_rho_g;
    val.gamma_s_over_rho_g = gamma_s_over_rho_g;
    val.gamma_phase_f = gamma_phase_f;
    val.gamma_phase_s = gamma_phase_s;
    val.w = w;
    val.phase_X_f = phase_X_f;
    val.phase_X_u = phase_X_u;
end

end