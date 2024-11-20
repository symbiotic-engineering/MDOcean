function [LCOE, P_var, P_matrix_elec, g, val] = simulation(X, p)	

X = max(X,1e-3); % sometimes optimizer will choose inputs that violate bounds, this is to prevent errors from negative numbers

%% Assemble inputs
in = p;
  
in.D_s   = X(1);     % inner diameter of float (m)
in.D_f   = X(2);     % outer diameter of float (m)
in.T_f_2 = X(3);     % draft of float (m)
in.h_s   = X(4);     % total height of spar (m)
in.F_max = X(5)*1e6; % max powertrain force (N)
in.B_p   = X(6)*1e6; % controller (powertrain) damping (Ns/m)
in.w_n   = X(7);     % controller (powertrain) natural frequency (rad/s)
in.M     = X(8);     % material (-)

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
    A_c, A_lat_sub, r_over_t, ...
    I, T, V_f_pct, V_s_pct, GM] = geometry(in.D_s, in.D_f, in.D_f_in, in.D_f_b, in.T_f_1, in.T_f_2, in.h_f, in.h_s, ...
                                            in.t_ft, in.t_fr, in.t_fc, in.t_fb, in.t_sr, ...
                                            in.t_dt, in.D_d, in.D_dt, in.theta_dt, in.T_s, in.h_d, in.t_d_max,...
                                            in.M, in.rho_m, in.rho_w, in.m_scale);

m_f_tot = max(m_f_tot,1e-3); % zero out negative mass produced by infeasible inputs

[F_heave_max, F_surge_max, F_ptrain_max, ...
	    P_var, P_avg_elec, P_matrix_elec, ...
        X_constraints] = dynamics(in, m_f_tot, m_s_tot, V_d, T);

[FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
                                    F_heave_max, F_surge_max,...
                                    in.M, in.h_s, in.T_s, in.rho_w, in.g, in.sigma_y, A_c, ...
                                    A_lat_sub, r_over_t, I, in.E);

LCOE = econ(m_m, in.M, in.cost_m, in.N_WEC, P_avg_elec, in.FCR, in.cost_perN, in.F_max, in.eff_array);

%% Assemble constraints g(x) >= 0
num_g = 15+numel(p.JPD);
g = zeros(1,num_g);
g(1) = V_f_pct;                         % prevent float too heavy
g(2) = 1 - V_f_pct;                     % prevent float too light
g(3) = V_s_pct;                         % prevent spar too heavy
g(4) = 1 - V_s_pct;                     % prevent spar too light
g(5) = GM;                              % pitch stability of float-spar system
g(6) = FOS1Y / p.FOS_min - 1;           % float survives max force
g(7) = FOS2Y / p.FOS_min - 1;           % spar survives max force
g(8) = FOS3Y / p.FOS_min - 1;           % damping plate survives max force
g(9) = FOS_buckling / p.FOS_min - 1;    % spar survives max force in buckling
g(10) = P_avg_elec;                     % positive power
%1 + min(Kp_over_Ks,[],'all');   % spar heave stability (positive effective stiffness)
g(11) = p.LCOE_max/LCOE - 1;            % prevent more expensive than threshold
g(12) = F_ptrain_max/in.F_max - 1;      % prevent irrelevant max force -
                                        % this constraint should always be active
                                        % and is only required when p.cost_perN = 0.
g(13) = X_constraints(1);               % prevent float rising above top of spar
g(14) = X_constraints(2);               % prevent float going below bottom of spar
g(15) = X_constraints(3);               % float amplitude obeys linear theory
g(16:end) = X_constraints(4:end);       % prevent rising out of water/slamming

criteria = all(~isinf(g)) && all(~isnan(g)) && all(isreal(g));
%assert( criteria )
if ~criteria
    warning('Inf, NaN, or imaginary constraint detected')
end

if nargout > 4 % if returning extra struct output for validation
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
     mass, CB_f,CG_f] = geometry(in.D_s, in.D_f, in.D_f_in, in.D_f_b, in.T_f_1, in.T_f_2,...
                                 in.h_f, in.h_s, in.t_ft, in.t_fr, in.t_fc, ...
                                 in.t_fb, in.t_sr, in.t_dt, in.D_d, in.D_dt, ...
                                 in.theta_dt, in.T_s, in.h_d, in.t_d_max,...
                                 in.M, in.rho_m, in.rho_w, in.m_scale);
    [~,capex,opex] = econ(m_m, in.M, in.cost_m, in.N_WEC, P_avg_elec, in.FCR, in.cost_perN, in.F_max, in.eff_array);
    [~, ~, ~, ~, ~, ~, ~, B_p,X,P_matrix_mech] = dynamics(in, m_f_tot, m_s_tot, V_d, T);
    val.mass_f  = mass(1);
    val.mass_vc = mass(2);
    val.mass_rp = mass(3);
    val.mass_tot = sum(mass);
    val.capex = capex;
    val.opex = opex;
    val.LCOE = LCOE;
    val.power_avg = P_avg_elec;
    val.power_max = max(P_matrix_elec,[],'all');
    val.force_heave = F_heave_max;
    val.force_ptrain = F_ptrain_max;
    val.FOS_b = FOS_buckling;
	val.c_v = P_var;
    val.B_p = B_p;
    val.X = X;
    val.P_mech = P_matrix_mech;
    val.CB_f = CB_f;
    val.CG_f = CG_f;
    val.vol_f = V_d(1);
    val.vol_s = V_d(2) + V_d(3);
end

end