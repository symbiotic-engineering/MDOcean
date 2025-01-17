function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
          	F_heave_storm, F_surge_storm, F_heave_op, F_surge_op, ...              % forces
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt,...        % bulk dimensions
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, h_d, A_dt, h_stiff, w_stiff,...% structural dimensions
            M, rho_w, g, sigma_y, sigma_e, E, nu)                                  % constants

    F_heave_peak = max(F_heave_storm,F_heave_op);
    F_surge_peak = max(F_surge_storm,F_surge_op);

    % inputs for both DLCs
    shared_inputs = {h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, h_d, A_dt, h_stiff, w_stiff, ...
            rho_w, g, E(M), nu(M)};

    % DLC 1: peak
    [FOS1Y, FOS2Y, FOS3Y, FOS_buckling] = structures_one_case(...
            F_heave_peak, F_surge_peak, sigma_y(M), shared_inputs{:});
    
    % DLC 2: endurance limit (long cycle fatigue)
    [FOS1Y(2), FOS2Y(2), FOS3Y(2), FOS_buckling(2)] = structures_one_case(...
            F_heave_op, F_surge_op, sigma_e(M), shared_inputs{:});
end

function [FOS1Y, FOS2Y, FOS3Y, FOS_spar_local] = structures_one_case(...
            F_heave, F_surge, sigma_max, ...
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt,...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, h_d, A_dt, h_stiff, w_stiff, ...
            rho_w, g, E, nu)
    
    depth = T_s; % max depth
    P_hydrostatic = rho_w * g * depth;

    %% Float plate stress
    [sigma_float_bot,...
     sigma_float_top] = float_plate_stress(D_f, D_f_in, F_heave, num_sections, ...
                                           t_bot, t_top, h_stiff, w_stiff, D_f_tu, nu);
    % ignoring top and side plates for now

    %% Spar Buckling calculation
    [FOS_spar,FOS_spar_local] = spar_combined_buckling(F_heave, E, I(2), h_s, D_s, A_c(2), t_s_r, ...
                                        P_hydrostatic, sigma_max, nu);

    %% Damping plate
    radial_stress_damping_plate = damping_plate_structures(F_heave, D_d, D_s,P_hydrostatic,t_d,A_dt,...
                                            theta_dt,L_dt,h_d,A_c,E);

    %% Factor of Safety (FOS) Calculations
    FOS1Y = sigma_max / sigma_float_bot;
    FOS2Y = FOS_spar;
    FOS3Y = sigma_max / radial_stress_damping_plate;

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

function sigma_vm = damping_plate_structures(F_heave, D_d, D_s,P_hydrostatic,t_d,A_dt,...
                                            theta_dt,L_dt,h_d,A_c,E)
    %% Stress calculations


    % calculate the deflection of the damping plate
    r = D_d/2;
    r_bending = r-(D_s/2);

    x = linspace(0,r_bending-0.01,1500);
    x1 = linspace(0,r-t_d,1500);

    [sigma_xx,shear,sigma_bending,...
     sigma_axial,M_x,I_x,F_water,...
     F_support,y] = damping_plate_func(E,D_d,A_dt,theta_dt,L_dt,h_d,D_s,h_d - t_d,...
                                          A_c(1),A_c(2),A_c(3),...
                                          P_hydrostatic,P_hydrostatic,P_hydrostatic,...
                                          x,x1,F_heave);
    [y_tip_Roark,l1,l2] = Roark_func(r-0.01,D_s/2,D_d/2,F_support*sin(theta_dt),...
                                    t_d,E,F_water/r_bending,D_s,0.26);

    %plots for debugging
    figure
    hold on
    %plot(x,M_x,'DisplayName','M')
    plot(x,I_x,'DisplayName','I')
    %plot(x, sigma_xx,'DisplayName','bending')
    %plot(x, sigma_axial,'DisplayName','axial')
    %plot(x,y)
    legend
    hold off

    %sigma_surge = F_surge ./ A_lat_sub;
    %sigma_rr = P_hydrostatic + sigma_surge;     % radial compression
    %sigma_tt = 0;%P_hydrostatic .* r_over_t;       % hoop stress
    sigma_zz = F_heave ./ A_c(3);                  % axial compression
    %sigma_rt = sigma_surge;                     % shear

    sigma_rr = max(abs(sigma_xx));
    sigma_tt = 0;
    sigma_zz = sigma_zz + P_hydrostatic;
    sigma_rt = max(abs(shear));
    sigma_tz = 0;
    sigma_zr = 0;

    
    % uncomment for debugging
    %     sigma = zeros(3,3,3);
    %     for j=1:3
    %     sigma(:,:,j) = [sigma_rr(j) sigma_rt(j) sigma_zr(j);
    %                     sigma_rt(j) sigma_tt(j) sigma_tz(j);
    %                     sigma_zr(j) sigma_tz(j) sigma_zz(j)];
    %     end

    % assume ductile material for now - need to use mohr's circle for concrete
    sigma_vm = von_mises(sigma_rr, sigma_tt, sigma_zz, sigma_rt, sigma_tz, sigma_zr);

end

function s_vm = von_mises(s_11, s_22, s_33, s_12, s_23, s_31)

    principal_term = 1/2 * ( (s_11 - s_22).^2 + (s_22 - s_33).^2 + (s_33 - s_11).^2 );
    shear_term = 3 * (s_12.^2 + s_23.^2 + s_31.^2);
    
    s_vm = sqrt( principal_term + shear_term );

end

