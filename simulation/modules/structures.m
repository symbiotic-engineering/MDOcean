function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
          	F_heave, F_surge, M, h_s, T_s, rho_w, g, ...
            sigma_y, A_c, A_lat_sub, r_over_t, I, E)

    %% Stress calculations
    depth = [0 T_s T_s]; % max depth
    
    P_hydrostatic = rho_w * g * depth;
    sigma_surge = F_surge ./ A_lat_sub;
    
    sigma_rr = P_hydrostatic + sigma_surge;     % radial compression
    sigma_tt = P_hydrostatic .* r_over_t;       % hoop stress
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
    
    %% Buckling calculation
    K = 2; % fixed-free - top is fixed by float angular stiffness, bottom is free
    L = h_s;
    F_buckling = pi^2 * E(M) * I(2) / (K*L)^2;
    
    %% Factor of Safety (FOS) Calculations
    FOS_yield = sigma_y(M) ./ sigma_vm;
    FOS1Y = FOS_yield(1);
    FOS2Y = FOS_yield(2);
    FOS3Y = FOS_yield(3);
    FOS_buckling = F_buckling ./ F_heave;

end 

function s_vm = von_mises(s_11, s_22, s_33, s_12, s_23, s_31)

    principal_term = 1/2 * ( (s_11 - s_22).^2 + (s_22 - s_33).^2 + (s_33 - s_11).^2 );
    shear_term = 3 * (s_12.^2 + s_23.^2 + s_31.^2);
    
    s_vm = sqrt( principal_term + shear_term );

end

