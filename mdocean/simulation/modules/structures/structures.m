function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
          	F_heave_storm, F_surge_storm, F_heave_op, F_surge_op, ...                       % forces
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt,D_d_tu,...% bulk dimensions
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, t_d_tu, h_d, A_dt, ...             % structural dimensions
            h_stiff_f, w_stiff_f, h_stiff_d, w_stiff_d, ...                                 % stiffener dimensions
            M, rho_w, g, sigma_y, sigma_e, E, nu, ...                                       % constants
            num_terms_plate, radial_mesh_plate, num_stiff_d)                                % plate hyperparameters

    F_heave_peak = max(F_heave_storm,F_heave_op);
    F_surge_peak = max(F_surge_storm,F_surge_op);

    % inputs for both DLCs
    shared_inputs = {h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt, D_d_tu,...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, t_d_tu, h_d, A_dt, ...
            h_stiff_f, w_stiff_f, h_stiff_d, w_stiff_d,...
            rho_w, g, E(M), nu(M), num_terms_plate, radial_mesh_plate, num_stiff_d};

    % DLC 1: peak
    sigma_buckle = sigma_y(M);     % fixme: to find ultimate, need to implement the ABS buckling formulas.
    sigma_u = sqrt(sigma_y(M) * sigma_buckle);
    [FOS1Y, FOS2Y, FOS3Y, FOS_buckling] = structures_one_case(...
            F_heave_peak, F_surge_peak, sigma_u, shared_inputs{:});
    
    % DLC 2: endurance limit (long cycle fatigue)
    [FOS1Y(2), FOS2Y(2), FOS3Y(2), FOS_buckling(2)] = structures_one_case(...
            F_heave_op, F_surge_op, sigma_e(M), shared_inputs{:});

end

function [FOS1Y, FOS2Y, FOS3Y, FOS_spar_local] = structures_one_case(...
            F_heave, F_surge, sigma_max, ...
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt, D_d_tu,...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, t_d_tu, h_d, A_dt, ...
            h_stiff_f, w_stiff_f, h_stiff_d, w_stiff_d,...
            rho_w, g, E, nu, num_terms_plate, radial_mesh_plate, num_stiff_d)
    
    depth = T_s; % max depth
    P_hydrostatic = rho_w * g * depth;

    %% Float plate stress
    [sigma_float_bot,...
     sigma_float_top] = float_plate_stress(D_f, D_f_in, F_heave, num_sections, ...
                                           t_bot, t_top, h_stiff_f, w_stiff_f, D_f_tu, nu);
    % ignoring top and side plates for now

    %% Spar Buckling calculation
    [FOS_spar,FOS_spar_local] = spar_combined_buckling(F_heave, E, I(2), h_s, D_s, A_c(2), t_s_r, ...
                                        P_hydrostatic, sigma_max, nu);

    %% Damping plate
    plot_on = false;
    sigma_damping_plate = damping_plate_structures(F_heave, D_d, D_s,P_hydrostatic,t_d,A_dt,...
                                            theta_dt,L_dt,h_d,A_c,E,nu, h_stiff_d,w_stiff_d, ...
                                            D_d_tu, t_d_tu, num_terms_plate, radial_mesh_plate, ...
                                            num_stiff_d, plot_on);

    %% Factor of Safety (FOS) Calculations
    FOS1Y = sigma_max / sigma_float_bot;
    FOS2Y = FOS_spar;
    FOS3Y = sigma_max / sigma_damping_plate;

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

function [sigma_vm,delta_max] = damping_plate_structures(F_heave, D_d, D_s,P_hydrostatic,t_d,A_dt,...
                                            theta_dt,L_dt,h_d,A_c,E,nu, h_stiff,width_stiff,...
                                            D_d_tu, t_d_tu, N, radial_mesh_plate, ...
                                            num_stiffeners, plot_on)
    
    a = D_d/2;
    b = D_s/2;

    rho = linspace(b/a, 1, radial_mesh_plate); % evaluate at these radial points

    % stiffness ratio and dimensional section properties
    [plate_tube_K_ratio, h_eq_vec, y_max_vec, D_eq] = get_plate_tube_stiffness(...
                                                    h_stiff, num_stiffeners, ...
                                                    rho, a, t_d, width_stiff, E, nu, ...
                                                    D_d_tu, t_d_tu, ...
                                                    L_dt, D_d, D_s);
    
    % nondimensional moment, displacement, and tube force ratio
    [Mr_over_F_heave, ...
     Mt_over_F_heave, ...
     Q_nondim, ...
     delta_over_a2Fheave_2piDeq, ...
     F_tube_over_F_heave]     = combined_plate_nondim(rho, b, a, nu, N, ...
                                                      plate_tube_K_ratio);

    % dimensional forces/moments
    F_tube = F_tube_over_F_heave      * F_heave;
    Mr     = Mr_over_F_heave / (2*pi) * F_heave;
    Mt     = Mt_over_F_heave / (2*pi) * F_heave;
    Q      = Q_nondim * a    / (2*pi) * F_heave;

    % dimensional displacement
    delta_total = a^2 * F_heave / (2*pi * D_eq) * delta_over_a2Fheave_2piDeq;
    delta_max = max(abs(delta_total));

    % dimensional stress
    sigma_rr_vec = get_plate_stress(Mr, y_max_vec, h_eq_vec);
    sigma_tt_vec = get_plate_stress(Mt, t_d/2,     t_d); % stiffeners don't affect section properties in theta 
    sigma_rz_vec = Q / t_d;

    sigma_zz = 0; sigma_rt = 0; sigma_tz = 0;
    sigma_vm_vec = von_mises(sigma_rr_vec, sigma_tt_vec, sigma_zz, ...
                             sigma_rt,     sigma_tz,     sigma_rz_vec);

    % maximum stress at all radii
    sigma_vm  = max(abs(sigma_vm_vec));

    if plot_on
        plot_damping_plate(a, D_eq, F_heave, F_tube, r, ...
            delta_plate_dis_nondim_vec, delta_plate_con_nondim_vec, ...
            Mr_con_nondim, Mr_dis_nondim_vec, Mr, ...
            y_max_vec, h_eq_vec, sigma_rr_vec)
    end
end

function [Mr_over_F_heave, ...
          Mt_over_F_heave, ...
          Q_nondim, ...
          delta_over_a2Fheave_2piDeq, ...
          F_tube_over_F_heave]     = combined_plate_nondim(rho, b, a, nu, N, plate_tube_K_ratio)

    % nondimensional annular plate solutions
    theta = [0 pi/2 pi 3*pi/2]; % evaluate at these angles - simulating a 
    % single concentrated load applied at theta=0 evaluated at these 4 angles 
    % and summing is the same as simulating 4 concentrated loads applied at
    % these angles evaluated at theta=0, due to superposition.

    % distrubuted load
    q = 1; % actually q = F_heave / (pi * (a^2-b^2)) but q=1 makes it nondimensional
    [delta_plate_dis_nondim_vec, ... % defined as delta_dis * 2pi * D / (F_heave * a^2) evaluated at all r's
     Mr_dis_nondim_vec, ...          % defined as Mr_dis * 2pi / F_heave evaluated at all r's
     Mt_dis_nondim_vec, ...          % defined as Mt_dis * 2pi / F_heave evaluated at all r's
     Q_dis_nondim_vec] ...           % defined as Q_dis * 2pi / (F_heave * a) evaluated at all r's
                        = distributed_plate_nondim(a,b,q,nu,rho);

    % concentrated load
    [delta_plate_con_nondim_vec, ... % defined as delta_con * 2pi * D / (F_tube * a^2) for a single F_tube evaluated at 4 thetas and all r's
     Mr_con_nondim_vec] ...          % defined as    Mr_con * 2pi / (F_tube)           for a single F_tube evaluated at 4 thetas and all r's
                        = concentrated_plate_nondim(b/a,nu,theta,rho,N);
    Mt_con_nondim = 0; % solutions aren't available for theta bending moment and shear for concentrated load, so assume zero
    Q_con_nondim = 0;

    % use deflection at outer edge for compatibility
    delta_plate_dis_nondim = delta_plate_dis_nondim_vec(end);
    delta_plate_con_nondim = delta_plate_con_nondim_vec(end,:);
    
    % sum the four angles for concentrated solution
    delta_plate_con_nondim = sum(delta_plate_con_nondim); % total for four equally-spaced F_tubes
    Mr_con_nondim = sum(Mr_con_nondim_vec,2);             % total for four equally-spaced F_tubes

    % compatibility using tube stiffness (derived from equating plate and tube deflection)
    % see notebook 7 p70 1/16/25 and notebook 9 p1 11/20/25
    F_tube_over_F_heave = delta_plate_dis_nondim / (plate_tube_K_ratio - delta_plate_con_nondim);

    % nondim moment superposition
    Mr_con_over_F_heave = Mr_con_nondim.' * F_tube_over_F_heave; % defined as Mr_con * 2 * pi / F_heave for four equally-spaced F_tubes
    Mr_dis_over_F_heave = Mr_dis_nondim_vec;                     % defined as Mr_dis * 2 * pi / F_heave
    Mr_over_F_heave = Mr_con_over_F_heave + Mr_dis_over_F_heave; % defined as Mr * 2 * pi / F_heave for combined loading (dist load F_heave and four equally-spaced point loads F_tube) 

    Mt_con_over_F_heave = Mt_con_nondim.' * F_tube_over_F_heave;
    Mt_dis_over_F_heave = Mt_dis_nondim_vec;
    Mt_over_F_heave = Mt_con_over_F_heave + Mt_dis_over_F_heave;

    % nondim shear force superposition
    Q_nondim = Q_con_nondim.' * F_tube_over_F_heave + Q_dis_nondim_vec;

    % nondim displacement superposition
    delta_over_a2Fheave_2piDeq = delta_plate_dis_nondim_vec ...
                                + F_tube_over_F_heave * sum(delta_plate_con_nondim_vec.',1); % defined as delta * 2pi * D / (F_heave * a^2) for combined loading
end

function [plate_tube_K_ratio, h_eq_vec, y_max_vec, D_eq] = get_plate_tube_stiffness(h_stiff,num_stiffeners, ...
                                                        rho, a, t_d, width_stiff, E, nu, ...
                                                        D_d_tu, t_d_tu, ...
                                                        L_dt, D_d, D_s)
    % stiffeners
    num_unique_stiffeners = length(h_stiff)/2;
    num_stiffener_repeats = num_stiffeners / num_unique_stiffeners;
    r = rho * a;
    circumf = 2*pi*r;
    width_plate = circumf / num_stiffener_repeats;

    % fixme should this be t_d or t_d*2? Technically should account for gap
    [h_eq_vec,y_max_vec] = get_stiffener_equivalent_properties(t_d, h_stiff, width_plate, width_stiff);
    h_eq = mean(h_eq_vec); % approximation to avoid dealing with radially varying equivalent thickness

    % equivalent D
    D_eq = E * h_eq^3 / (12*(1-nu^2));

    % tube bending stiffness - p72 notebook 1/17/25
    I_tube = pi/64 * ( D_d_tu^4 - (D_d_tu - 2*t_d_tu)^4 );
    K_tube = 6*E*I_tube / (L_dt^2 * (D_d - D_s));

    % notebook p1 11/20/25
    plate_tube_K_ratio = D_eq*2*pi/(a^2*K_tube);
end

function [] = plot_damping_plate(r, delta_total, ...
            delta_plate_dis_nondim_vec, delta_plate_con_nondim_vec, ...
            Mr_con_nondim, Mr_dis_nondim_vec, Mr, ...
            y_max_vec, h_eq_vec, sigma_r_vec)
    
    figure
    plot(r,Mr_con_nondim,'DisplayName','Mr con nondim')
    hold on
    plot(r,Mr_dis_nondim_vec,'DisplayName','Mr dis nondim')
    plot(r,Mr/max(abs(Mr)),'DisplayName','Mr normalized')
    legend
    xlabel('r')
    ylabel('Moment')
    improvePlot

    figure
    plot(r,delta_plate_dis_nondim_vec,'DisplayName','delta dis nondim')
    hold on
    plot(r,delta_plate_con_nondim_vec(:,1),'DisplayName','delta con nondim: theta=0')
    plot(r,sum(delta_plate_con_nondim_vec,2),'DisplayName','delta con nondim: sum all 4 theta')
    plot(r,delta_total/max(abs(delta_total)),'DisplayName','delta total normalized')
    legend
    xlabel('r')
    ylabel('Deflection')
    improvePlot

    figure
    plot(r,y_max_vec/max(y_max_vec),'DisplayName','y max normalized')
    hold on
    plot(r,h_eq_vec/max(h_eq_vec),'DisplayName','h eq normalized')
    plot(r,sigma_r_vec/sigma_r,'DisplayName','sigma r normalized')
    legend
    xlabel('r')
    ylabel('Normalized Quantity')
    improvePlot

end

function s_vm = von_mises(s_11, s_22, s_33, s_12, s_23, s_31)

    principal_term = 1/2 * ( (s_11 - s_22).^2 + (s_22 - s_33).^2 + (s_33 - s_11).^2 );
    shear_term = 3 * (s_12.^2 + s_23.^2 + s_31.^2);
    
    s_vm = sqrt( principal_term + shear_term );

end

