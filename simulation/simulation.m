function [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, D_d, P_matrix, g, val] = simulation(X, p)	

%% Assemble inputs
in = p;

in.D_f              = X(1);     % outer diameter of float (m)
D_s_over_D_f        = X(2);     % normalized diameter of spar column (-)
h_f_over_D_f        = X(3);     % normalized vertical thickness of float (-)
T_s_over_h_s        = X(4);     % normalized spar draft below waterline (-)
in.F_max            = X(5)*1e6; % max powertrain force (N)
in.B_p              = X(6)*1e6; % controller (powertrain) damping (Ns/m)
in.w_n              = X(7);     % controller (powertrain) natural frequency (rad/s)
in.M                = X(8);     % material (-)

% Variable ratios defined by design variables
in.D_s = D_s_over_D_f * in.D_f;
in.h_f = h_f_over_D_f * in.D_f;
% Geometric similarity to float submergence parameter
in.T_f = p.T_f_over_h_f * in.h_f;
% Geometric similarity to maintain constant damping ratio
% D_s sets D_d, T_s, h_d
D_d = p.D_d_over_D_s * in.D_s;
in.T_s = p.T_s_over_D_s * in.D_s;
in.h_d = p.h_d_over_D_s * in.D_s;
% Another ratio defined by design variable
in.h_s = 1/T_s_over_h_s * in.T_s;


%% Run modules
[V_d, m_m, m_f_tot, ...
    A_c, A_lat_sub, r_over_t, ...
    I, T, V_f_pct, V_s_pct, GM] = geometry(in.D_s, in.D_f, in.T_f, in.h_f, in.h_s, ...
                                            in.t_ft, in.t_fr, in.t_fc, in.t_fb, in.t_sr, ...
                                            D_d, in.T_s, in.h_d, ...
                                            in.M, in.rho_m, in.rho_w, in.m_scale);
B = [V_f_pct, V_s_pct]; % temporary to avoid changing output of simulation

[F_heave_max, F_surge_max, F_ptrain_max, ...
	    P_var, P_elec, P_matrix, h_s_extra] = dynamics(in, m_f_tot, V_d, T);

[FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
                                    F_heave_max, F_surge_max,...
                                    in.M, in.h_s, in.T_s, in.rho_w, in.g, in.sigma_y, A_c, ...
                                    A_lat_sub, r_over_t, I, in.E);

LCOE = econ(m_m, in.M, in.cost_m, in.N_WEC, P_elec, in.FCR, in.array_eff);

%% Assemble constraints g(x) >= 0
g = zeros(1,18);
g(1) = V_f_pct;                         % prevent float too heavy
g(2) = 1 - V_f_pct;                     % prevent float too light
g(3) = V_s_pct;                         % prevent spar too heavy
g(4) = 1 - V_s_pct;                     % prevent spar too light
g(5) = GM;                              % stability
g(6) = FOS1Y(1) - p.FOS_min;            % float survives hydro force
%g(7) = FOS1Y(2) - p.FOS_min;            % float survives powertrain force
g(8) = FOS2Y(1) - p.FOS_min;            % spar survives hydro force
%g(9) = FOS2Y(2) - p.FOS_min;            % spar survives powertrain force
g(10) = FOS_buckling(1) - p.FOS_min;    % spar survives hydro force in buckling
%g(11) = FOS_buckling(2) - p.FOS_min;    % spar survives powertrain force in buckling
g(12) = FOS3Y(1) - p.FOS_min;           % damping plate survives hydro force
%g(13) = FOS3Y(2) - p.FOS_min;           % damping plate survives powertrain force
g(14) = P_elec;                         % positive power
g(15) = D_d - p.D_d_min;                % damping plate diameter (spar natural freq)
g(16) = h_s_extra;                      % prevent float rising above top of spar
g(17) = p.LCOE_max/LCOE - 1;            % prevent more expensive than threshold
g(18) = F_ptrain_max/in.F_max - 1;      % prevent irrelevant max force

assert( all(~isinf(g)) && all(~isnan(g)) )

if nargout > 12 % if returning extra struct output for validation
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, mass] = geometry(in.D_s, in.D_f, in.T_f, ...
                                            in.h_f, in.h_s, in.t_ft, in.t_fr, ...
                                            in.t_fc, in.t_fb, in.t_sr, ...
                                            D_d, in.T_s, in.h_d, ...
                                            in.M, in.rho_m, in.rho_w, in.m_scale);
    [~,capex,opex] = econ(m_m, in.M, in.cost_m, in.N_WEC, P_elec, in.FCR, in.array_eff);
    val.mass_f  = mass(1);
    val.mass_vc = mass(2);
    val.mass_rp = mass(3);
    val.capex = capex;
    val.opex = opex;
    val.LCOE = LCOE;
    val.power_avg = P_elec;
    val.power_max = max(P_matrix,[],'all');
    val.force_heave = F_heave_max;
    val.FOS_b = FOS_buckling;
end

end