function [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, P_elec, P_matrix] = simulation(X, p)	
        
% capital X is design variables in vector format (necessary for optimization)	
% lowercase x is design variables in struct format (more readable)	
x = struct( 'D_sft',X(1),...            % outer diameter of float (m)
            'D_i',  X(2)*X(1),...       % inner diameter of float (m)
            't_f',  X(3)*X(1),...       % vertical thickness of float (m)
            'd_f',  X(4)*X(3)*X(1),...  % float draft below waterline (m)
            'F_max',X(5)*1e6,...        % max powertrain force (N)
            'D_int',X(6)*1e6,...        % internal damping of controller (Ns/m)
            'w_n',  X(7),...            % internal spring of controller (N/m)
            'M',    X(8) );             % material (-)

% merge design variables and parameters to single input struct
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
in = mergestructs(x,p);

[V_d, V_m, m_tot, m_float, h,...
    A_c, A_lat_sub, r_over_t, I, draft] = geometry(in.D_i, in.D_sft, in.t_sft, ...
                                            in.t_sf, in.t_sfb, in.t_vc, ...
                                            in.D_or, in.t_r, in.rho_m, in.M, in.rho_w);
        	
[F_hydro_heave, F_hydro_surge, ...
	F_ptrain, P_var, P_elec, P_matrix] = dynamics(x, p, m_float, V_d, draft);

[B,FOS1Y,FOS2Y,FOS3Y,FOS_buckling,GM] = structures(...
                                    F_hydro_heave, F_hydro_surge, F_ptrain,...
                                    in.M, h, in.rho_w, in.g, in.sigma_y, A_c, ...
                                    A_lat_sub, r_over_t, I, in.E);

LCOE = econ(m_tot, V_m, in.M, in.cost_m, in.N_WEC, P_elec, in.FCR);

end