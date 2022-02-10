function [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, P_matrix] = simulation(X, p)	
        
% capital X is design variables in vector format (necessary for optimization)	
% lowercase x is design variables in struct format (more readable)	
x = struct( 'D_sft',X(1),...        % outer diameter of float (m)
            'D_i',  X(2)*X(1),...   % inner diameter of float (m)
            'D_or', X(3)*X(1),...        % outer diameter of reaction plate (m)
            'M',    X(4),...        % material (-)
            'N_WEC',X(5),...        % number of WECs in array (-)
            'D_int',X(6)*1e6,...    % internal damping of controller (Ns/m)
            'w_n',  X(7) );         % internal spring of controller (N/m)
        
[V_d, V_m, m_tot, m_float, h, t_f,...
    A_c, A_lat_sub, r_over_t, I, draft] = geometry(x.D_i, x.D_sft, p.t_sft, ...
                                            p.t_sf, p.t_sfb, p.t_vc, ...
                                            x.D_or, p.t_r, p.rho_m, x.M);
        	
[F_hydro_heave, F_hydro_surge, ...
	F_ptrain, P_var, P_elec, P_matrix] = dynamics(x, p, m_float, V_d, draft);

[B,FOS1Y,FOS2Y,FOS3Y,FOS_buckling,GM] = structures(V_d, m_tot, ...
                                    F_hydro_heave, F_hydro_surge, F_ptrain,...
                                    x.M, h, p.rho_w, p.g, p.sigma_y, A_c, ...
                                    A_lat_sub, r_over_t, I, p.E, t_f, p.t_r);

LCOE = econ(m_tot, V_m, x.M, p.cost_m, x.N_WEC, P_elec, p.FCR);

end