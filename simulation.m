
function [LCOE, D_env, Lt, B, FOS, P_elec] = simulation(X, p)

% capital X is design variables in vector format (necessary for optimization)
% lowercase x is design variables in struct format (more readable)

x = struct( 'D_sft',X(1),...        % outer diameter of float (m)
            'D_i',  X(2)*X(1),...   % inner diameter of float (m)
            'L_sf', pi*X(1),...        % radial material thickness of float (m) 
            'D_or', X(3),...        % outer diameter of reaction plate (m)
            'M',    X(4),...        % material (-)
            'N_WEC',X(5),...        % number of WECs in array (-)
            'D_int',X(6));          % internal damping of controller (Ns/m)

[V_d, V_m, m_tot, m_float, h, ...
    t_f, A_c, A_lat_sub, r_over_t, I] = geometry(x.D_i, x.D_sft, p.t_sft, ...
                                        x.L_sf, p.t_sf, p.t_sfb, p.t_vc, ...
                                        x.D_or, p.t_r, p.rho_m, x.M);
        
[F_hydro_heave, F_hydro_surge, ...
        F_ptrain, D_env, P_elec] = dynamicSimulation(x, p, m_float, t_f);

[B,FOS] = structures(V_d, m_tot, F_hydro_heave, F_hydro_surge, F_ptrain, ...
        x.M, h, p.rho_w, p.g, p.sigma_y, A_c, A_lat_sub, r_over_t, I, p.E);

[LCOE, Lt] = econ(m_tot, V_m, x.M, p.cost_m, x.N_WEC, P_elec);

end