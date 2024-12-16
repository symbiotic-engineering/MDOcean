function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
          	F_heave_storm, F_surge_storm, F_heave_op, F_surge_op, ... % forces
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, ... % bulk dimensions
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, h_stiff, w_stiff, ... % structural dimensions
            M, rho_w, g, sigma_y, sigma_e, E, nu) % constants

    F_heave_peak = max(F_heave_storm,F_heave_op);
    F_surge_peak = max(F_surge_storm,F_surge_op);

    % peak
    [FOS1Y, FOS2Y, FOS3Y, FOS_buckling] = structures_one_case(...
            F_heave_peak, F_surge_peak, sigma_y(M), ...
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, ...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, h_stiff, w_stiff, ...
            M, rho_w, g, E, nu);
    
    % endurance limit (fatigue)
    [FOS1Y(2), FOS2Y(2), FOS3Y(2), FOS_buckling(2)] = structures_one_case(...
            F_heave_op, F_surge_op, sigma_e(M), ...
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, ...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, h_stiff, w_stiff, ...
            M, rho_w, g, E, nu);
end

function [FOS1Y, FOS2Y, FOS3Y, FOS_spar_local] = structures_one_case(...
            F_heave, F_surge, sigma_max, ...
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, ...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, h_stiff, w_stiff, ...
            M, rho_w, g, E, nu)

    %% Stress calculations
    depth = [0 T_s T_s]; % max depth
    
    P_hydrostatic = rho_w * g * depth;
    sigma_surge = F_surge ./ A_lat_sub;
    
    sigma_rr = P_hydrostatic + sigma_surge;     % radial compression
    sigma_tt = 0;%P_hydrostatic .* r_over_t;       % hoop stress
    sigma_zz = F_heave ./ A_c;                  % axial compression
    sigma_rt = sigma_surge;                     % shear
    sigma_tz = [0 0 0];
    sigma_zr = [0 0 0];
    
    % uncomment for debugging
    %     sigma = zeros(3,3,3);
    %     for j=1:3
    %     sigma(:,:,j) = [sigma_rr(j) sigma_rt(j) sigma_zr(j);
    %                     sigma_rt(j) sigma_tt(j) sigma_tz(j);
    %                     sigma_zr(j) sigma_tz(j) sigma_zz(j)];
    %     end
    
    % assume ductile material for now - need to use mohr's circle for concrete
    sigma_vm = von_mises(sigma_rr, sigma_tt, sigma_zz, sigma_rt, sigma_tz, sigma_zr);
    
    %% Float plate stress
    [sigma_float_bot,sigma_float_top] = float_plate_stress(D_f, D_f_in, F_heave, num_sections, t_bot, t_top, h_stiff, w_stiff, D_f_tu, nu);
    % ignoring top and side plates for now

    %% Buckling calculation
    [FOS_spar,FOS_spar_local] = spar_combined_buckling(F_heave, E(M), I(2), h_s, D_s, A_c(2), t_s_r, ...
                                        P_hydrostatic(2), sigma_max, nu(M));

    %% Factor of Safety (FOS) Calculations
    FOS1Y = sigma_max / sigma_float_bot;
    FOS2Y = FOS_spar;
    FOS3Y = sigma_vm(3);

end 

function [FOS_spar, FOS_spar_local] = spar_combined_buckling(F, E, I, L, D, A, t, q, sigma_0, nu)
    % euler buckling
    K = 2; % fixed-free - top is fixed by float angular stiffness, bottom is free
    F_buckling = pi^2 * E * I / (K*L)^2;
    sigma_EA = F_buckling / A;

    % hoop stress
    sigma_theta = q*D/(2*t);

    % local buckling of plate element: 2/9.7
    k_s = 1.33; % uniform compression, fixed-free, from table 3 on page 30
    s = D; % fixme - not sure if this is correct
    P_r = 0.6; % steel proportional linear elastic limit
    sigma_Ex = k_s * pi^2 * E / (12 * (1-nu^2)) * (t/s)^2;
    if sigma_Ex <= P_r * sigma_0
        sigma_Cx = sigma_Ex;
    else
        sigma_Cx = sigma_0 * (1 - P_r * (1-P_r) * sigma_0/sigma_Ex);
    end

    % stress at failure: 2/7.3
    compact = true; % assumption, eventually this should be a calculation
    if compact
        sigma_F = sigma_0; % yield
    else
        sigma_F = sigma_Cx; % local buckling of plate element
    end

    % whether the failure mode is pure euler buckling or combined loading
    sigma_EA_thresh = P_r * sigma_F * (1 - sigma_theta / sigma_F);
    pure_buckling = sigma_EA <= sigma_EA_thresh; 
    if pure_buckling
        sigma_C_A_theta = sigma_EA;
    else
        zeta = 1 - P_r*(1 - P_r)* sigma_F/sigma_EA - sigma_theta/sigma_F;
        omega = 0.5 * sigma_theta/sigma_F * (1 - .5 * sigma_theta/sigma_F);
        Lambda = 1/2 * (zeta + sqrt(zeta^2 + 4*omega));
        sigma_C_A_theta = sigma_F * Lambda;
    end
    
    sigma_ac = F ./ A + q;
    FOS_spar = sigma_C_A_theta / sigma_ac;

    % local buckling of tube: final part of 2/7.3, which combines 2/9.1
    % and 2/9.5
    FOS_spar_local = 3; % fixme: need to implement
end

function s_vm = von_mises(s_11, s_22, s_33, s_12, s_23, s_31)

    principal_term = 1/2 * ( (s_11 - s_22).^2 + (s_22 - s_33).^2 + (s_33 - s_11).^2 );
    shear_term = 3 * (s_12.^2 + s_23.^2 + s_31.^2);
    
    s_vm = sqrt( principal_term + shear_term );

end

