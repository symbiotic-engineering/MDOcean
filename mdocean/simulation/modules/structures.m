function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
          	F_heave, F_surge, M, h_s, T_s, rho_w, g, ...
            sigma_y, A_c, A_lat_sub, r_over_t, I, E, D_d, A_dt, theta_dt, L_dt, h_d, D_s, t_d)

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

    F_water = P_hydrostatic(2) * A_c(3)/2;
    r = D_d/2;
    I_d = (1/12) * r * h_d^3;
    F_support = 3*F_water*r^3*A_dt/(8*sin(theta_dt)*(3*L_dt*I_d - (r^3*A_dt)));

    x1 = linspace(0,r-t_d,1500);
    A = 2*sqrt(r^2-x1.^2) * h_d - (2*sqrt(r^2-x1.^2) * t_d^3);
    sigma_axial = -F_support * cos(theta_dt) ./ A;
    v = abs(-F_water - (F_support*(sin(theta_dt))) + (F_water*x1/r));
    shear = max(v./A);

    %count = 1;
    % for x = 0:0.01:r-0.01 %avoid infinite moment of inertia
    %     i = count;
    %     M_x(i) = -(F_water+F_support*sin(theta_dt))*x + (F_water*x^2/(2*r)) + ...
    %     (F_water*r/2) + (F_support*r*sin(theta_dt));
    %     I_x(i) = (1/12) * 2*sqrt(r^2-x^2) * h_d^3;
    %     count = count+1;
    % end

    r_bending = r-D_s/2;
    x = linspace(0,r_bending-0.01,1500);
    M_x = -(F_water+F_support*sin(theta_dt))*x + (F_water.*x.^2/(2*r_bending)) + ...
        (F_water*r_bending/2) + (F_support*r_bending*sin(theta_dt));
    I_x = (1/12) * 2*sqrt(r_bending^2-x.^2) * h_d^3 - ((1/12) * 2*sqrt(r_bending^2-x.^2) * t_d^3);

    sigma_bending = abs(M_x .* h_d./(2.*I_x));
    sigma_xx = max(sigma_axial) + max(sigma_bending);

    % uncomment for debugging
    %     sigma = zeros(3,3,3);
    %     for j=1:3
    %     sigma(:,:,j) = [sigma_rr(j) sigma_rt(j) sigma_zr(j);
    %                     sigma_rt(j) sigma_tt(j) sigma_tz(j);
    %                     sigma_zr(j) sigma_tz(j) sigma_zz(j)];
    %     end
    
    % assume ductile material for now - need to use mohr's circle for concrete
    sigma_rr(3) = sigma_xx;
    sigma_tt(3) = 0;
    sigma_zz(3) = sigma_zz(3)+P_hydrostatic(3);
    sigma_rt(3) = shear;

    sigma_vm = von_mises(sigma_rr, sigma_tt, sigma_zz, sigma_rt, sigma_tz, sigma_zr);
    %sigma_vm = von_mises(sigma_xx, [0 0 0], sigma_zz+P_hydrostatic, shear, [0 0 0], [0 0 0]);
    
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