function [LCOE, P_var, P_matrix_elec, g, val] = simulation(X, p)	

X = max(X,1e-3); % sometimes optimizer will choose inputs that violate bounds, this is to prevent errors from negative numbers

%% Assemble inputs
in = p;

in.D_s   = X(1);      % inner diameter of float (m)
in.D_f   = X(2);      % outer diameter of float (m)
in.T_f_2 = X(3);      % draft of float (m)
in.h_s   = X(4);      % total height of spar (m)
in.F_max = X(5)*1e6;  % max powertrain force (N)
in.B_p   = X(6)*1e6;  % controller (powertrain) damping (Ns/m)
in.h_fs_clear = X(7); % vertical clearance between float tubes and spar when at rest (m)
in.M     = X(8);      % material (-)

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
in.D_f_in = p.D_f_in_over_D_s * in.D_s;

%% Run modules
[V_d, m_m, m_f_tot, m_s_tot, ...
 A_c, A_lat_sub, ...
 I, T, V_f_pct, V_s_pct, GM] = geometry(in.D_s, in.D_f, in.D_f_in, in.D_f_b, ...
                                        in.T_f_1, in.T_f_2, in.h_f, in.h_s, ...
                                        in.h_fs_clear, in.D_f_tu, in.t_f_t, ...
                                        in.t_f_r, in.t_f_c, in.t_f_b, in.t_f_tu, ...
                                        in.t_s_r, in.t_d_tu, in.D_d, in.D_d_tu, ...
                                        in.theta_d_tu, in.T_s, in.h_d, in.t_d_max,...
                                        in.M, in.rho_m, in.rho_w, in.m_scale);

m_f_tot = max(m_f_tot,1e-3); % zero out negative mass produced by infeasible inputs

[F_heave_storm, F_surge_storm, ...
 F_heave_op, F_surge_op, F_ptrain_max, ...
 P_var, P_avg_elec, P_matrix_elec, ...
                    X_constraints] = dynamics(in, m_f_tot, m_s_tot, V_d, T);

[FOS_float,FOS_spar,FOS_damping_plate,...
    FOS_spar_local] = structures(F_heave_storm, F_surge_storm, F_heave_op, F_surge_op, ...
                                 in.M, in.h_s, in.T_s, in.rho_w, in.g, in.sigma_y, in.sigma_e, ...
                                 A_c, A_lat_sub, in.D_s, in.t_s_r, I, in.E, in.nu);

LCOE = econ(m_m, in.M, in.cost_m, in.N_WEC, P_avg_elec, in.FCR, in.cost_perN, in.F_max, in.eff_array);

%% Assemble constraints g(x) >= 0
num_g = 19+numel(p.JPD);
g = zeros(1,num_g);
g(1) = V_f_pct;                         % prevent float too heavy
g(2) = 1 - V_f_pct;                     % prevent float too light
g(3) = V_s_pct;                         % prevent spar too heavy
g(4) = 1 - V_s_pct;                     % prevent spar too light
g(5) = GM;                              % pitch stability of float-spar system
g(6) = FOS_float(1) / p.FOS_min - 1;    % float survives max force
g(7) = FOS_float(2) / p.FOS_min - 1;    % float survives fatigue
g(8) = FOS_spar(1) / p.FOS_min - 1;           % spar survives max force
g(9) = FOS_spar(2) / p.FOS_min - 1;           % spar survives fatigue
g(10) = FOS_damping_plate(1) / p.FOS_min - 1; % damping plate survives max force
g(11) = FOS_damping_plate(2) / p.FOS_min - 1; % damping plate survives fatigue
g(12) = FOS_spar_local(1) / p.FOS_min - 1;    % spar survives max force in local buckling
g(13) = FOS_spar_local(2) / p.FOS_min - 1;    % spar survives fatigue in local buckling
g(14) = P_avg_elec;                     % positive power
%1 + min(Kp_over_Ks,[],'all');   % spar heave stability (positive effective stiffness)
g(15) = p.LCOE_max/LCOE - 1;            % prevent more expensive than threshold
g(16) = F_ptrain_max/in.F_max - 1;      % prevent irrelevant max force -
                                        % this constraint should always be active
                                        % and is only required when p.cost_perN = 0.
g(17) = X_constraints(1);               % prevent float rising above top of spar
g(18) = X_constraints(2);               % prevent float going below bottom of spar
g(19) = X_constraints(3);               % float amplitude obeys linear theory
g(20:end) = X_constraints(4:end);       % prevent rising out of water/slamming
% fixme: add another X_constraint that h_fs_clear is greater than X_u

criteria = all(~isinf(g)) && all(~isnan(g)) && all(isreal(g));
%assert( criteria )
if ~criteria
    warning('Inf, NaN, or imaginary constraint detected')
end

if nargout > 4 % if returning extra struct output for validation
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
     mass, CB_f,CG_f] = geometry(in.D_s, in.D_f, in.D_f_in, in.D_f_b, ...
                                 in.T_f_1, in.T_f_2, in.h_f, in.h_s, ...
                                 in.h_fs_clear, in.D_f_tu, in.t_f_t, ...
                                 in.t_f_r, in.t_f_c, in.t_f_b, in.t_f_tu, ...
                                 in.t_s_r, in.t_d_tu, in.D_d, in.D_d_tu, ...
                                 in.theta_d_tu, in.T_s, in.h_d, in.t_d_max,...
                                 in.M, in.rho_m, in.rho_w, in.m_scale);
    [~,capex,opex] = econ(m_m, in.M, in.cost_m, in.N_WEC, P_avg_elec, in.FCR, in.cost_perN, in.F_max, in.eff_array);
    [~, ~, ~, ~, ~, ~, ~, B_p,X_u,X_f,X_s,P_matrix_mech] = dynamics(in, m_f_tot, m_s_tot, V_d, T);
    val.mass_f  = mass(1);
    val.mass_vc = mass(2);
    val.mass_rp = mass(3);
    val.mass_tot = sum(mass);
    val.capex = capex;
    val.opex = opex;
    val.LCOE = LCOE;
    val.power_avg = P_avg_elec;
    val.power_max = max(P_matrix_elec,[],'all');
    val.force_heave = F_heave_storm;
    val.force_ptrain = F_ptrain_max;
    val.FOS_spar = FOS_spar(1);
	val.c_v = P_var;
    val.B_p = B_p;
    val.X_u = X_u;
    val.X_f = X_f;
    val.X_s = X_s;
    val.P_mech = P_matrix_mech;
    val.CB_f = CB_f;
    val.CG_f = CG_f;
    val.vol_f = V_d(1);
    val.vol_s = V_d(2) + V_d(3);
end

end