function [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, D_d, P_matrix, g] = simulation(X, p)	

%% Assemble inputs
in = p;

in.D_f              = X(1);     % outer diameter of float (m)
D_s_over_D_f        = X(2);     % normalized diameter of spar column (-)
h_f_over_D_f        = X(3);     % normalized vertical thickness of float (-)
T_s_over_h_s        = X(4);     % normalized spar draft below waterline (-)
in.F_max            = X(5)*1e6; % max powertrain force (N)
in.D_int            = X(6)*1e6; % internal damping of controller (Ns/m)
in.w_n              = X(7);     % internal spring of controller (N/m)
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
                                            in.M, in.rho_m, in.rho_w);
B = [V_f_pct, V_s_pct]; % temporary to avoid changing output of simulation
        	
[F_hydro_heave, F_hydro_surge, F_ptrain, ...
	    P_var, P_elec, P_matrix, h_s_extra] = dynamics(in, m_f_tot, V_d, T);

[FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
                                    F_hydro_heave, F_hydro_surge, F_ptrain,...
                                    in.M, in.h_s, in.rho_w, in.g, in.sigma_y, A_c, ...
                                    A_lat_sub, r_over_t, I, in.E);

LCOE = econ(m_m, in.M, in.cost_m, in.N_WEC, P_elec, in.FCR);

%% Assemble constraints g(x) >= 0
g = zeros(1,14);
g(1) = V_f_pct;                         % prevent float too heavy
g(2) = 1 - V_f_pct;                     % prevent float too light
g(3) = V_s_pct;                         % prevent spar too heavy
g(4) = 1 - V_s_pct;                     % prevent spar too light
g(5) = GM;                              % stability
g(6) = FOS1Y(1) - p.FOS_min;            % float survives hydro force
g(7) = FOS1Y(2) - p.FOS_min;            % float survives powertrain force
g(8) = FOS2Y(1) - p.FOS_min;            % spar survives hydro force
g(9) = FOS2Y(2) - p.FOS_min;            % spar survives powertrain force
g(10) = FOS_buckling(1) - p.FOS_min;    % spar survives hydro force in buckling
g(11) = FOS_buckling(2) - p.FOS_min;    % spar survives powertrain force in buckling
g(12) = FOS3Y(1) - p.FOS_min;           % damping plate survives hydro force
g(13) = FOS3Y(2) - p.FOS_min;           % damping plate survives powertrain force
g(14) = P_elec;                         % positive power
g(15) = D_d - p.D_d_min;                % damping plate diameter (spar natural freq)
g(16) = h_s_extra;                      % prevent float rising above top of spar
g(17) = p.LCOE_max/LCOE - 1;            % prevent more expensive than nominal

assert( all(~isinf(g)) && all(~isnan(g)) )

end