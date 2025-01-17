function [FOS1Y,FOS2Y,FOS3Y,FOS_buckling] = structures(...
          	F_heave_storm, F_surge_storm, F_heave_op, F_surge_op, ...              % forces
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt,D_d_tu,...        % bulk dimensions
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, t_d_tu, h_d, A_dt, h_stiff, w_stiff,...% structural dimensions
            M, rho_w, g, sigma_y, sigma_e, E, nu)                                  % constants

    F_heave_peak = max(F_heave_storm,F_heave_op);
    F_surge_peak = max(F_surge_storm,F_surge_op);

    % inputs for both DLCs
    shared_inputs = {h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt, D_d_tu,...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, t_d_tu, h_d, A_dt, h_stiff, w_stiff, ...
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
            h_s, T_s, D_s, D_f, D_f_in, num_sections, D_f_tu, D_d, L_dt, theta_dt, D_d_tu,...
            t_s_r, I, A_c, A_lat_sub, t_bot, t_top, t_d, t_d_tu, h_d, A_dt, h_stiff, w_stiff, ...
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
                                            theta_dt,L_dt,h_d,A_c,E,nu, h_stiff,w_stiff, D_d_tu, t_d_tu);

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
                                            theta_dt,L_dt,h_d,A_c,E,nu, h_stiff,width_stiff,...
                                            D_d_tu, t_d_tu)
    
    a = D_d/2;
    b = D_s/2;
    N = 100;

     % nondimensional annular plate solutions
    rho = linspace(b/a,1,20); % evaluate at these radial points
    theta = [0 pi/2 pi 3*pi/2];
    [delta_plate_dis_nondim_vec, Mr_dis_nondim_vec, ...
                                 Mt_dis_nondim_vec] = distributed_plate_nondim(a,b,F_heave,nu,rho);
    [delta_plate_con_nondim_vec, Mr_con_nondim_vec] = concentrated_plate_nondim(b/a,nu,theta,rho,N);

    % max deflection at outer edge, 
    delta_plate_dis_nondim = abs(delta_plate_dis_nondim_vec(end));
    delta_plate_con_nondim = abs(delta_plate_con_nondim_vec(end,:));
    % max moment at inner edge
    %Mr_dis_nondim = abs(Mr_dis_nondim_vec(1));
    %Mr_con_nondim = abs(Mr_con_nondim_vec(1,:));
    
    % sum the four angles for concentrated solution
    delta_plate_con_nondim = sum(delta_plate_con_nondim);
    Mr_con_nondim = sum(Mr_con_nondim_vec,2);

    % stiffeners
    num_stiffeners = 24;
    r = rho * a;
    circumf = 2*pi*r;
    width_plate = circumf / num_stiffeners;

    % fixme should this be t_d or t_d*2? Technically should account for gap
    [h_eq_vec,y_max_vec] = get_stiffener_equivalent_properties(t_d, h_stiff, width_plate, width_stiff);
    h_eq = mean(h_eq_vec); % approximation to avoid dealing with radially varying equivalent thickness

    % equivalent D
    D_eq = E * h_eq^3 / (12*(1-nu^2));

    % tube bending stiffness - p72 notebook 1/17/25
    I_tube = pi/64 * ( D_d_tu^4 - (D_d_tu - 2*t_d_tu)^4 );
    K_tube = 6*E*I_tube / (L_dt^2 * (D_d - D_s));

    % compatbility (derived from equating plate and tube deflection)
    F_tube = F_heave * delta_plate_dis_nondim / (D_eq/(a^2*K_tube) - delta_plate_con_nondim);

    % moment superposition
    Mr_con = Mr_con_nondim.' * F_tube;
    Mr_dis = Mr_dis_nondim_vec * F_heave;
    Mr = Mr_con + Mr_dis;

    % stress
    sigma_r_vec = get_plate_stress(Mr, y_max_vec, h_eq_vec);
    sigma_r  = max(abs(sigma_r_vec));

    sigma_vm  = sigma_r; % ignore Mt for now
end

function [w_nondim,Mr_nondim,Mt_nondim] = distributed_plate_nondim(a,b,F_heave,nu,rho)
    A = pi*a^2;
    q = F_heave/A;
    P = F_heave/4;
    v = nu;
    r = rho*a;
    r0 = b;

    C2 = 1/4*(1-(b/a)^2*(1+2*log(a/b)));
    C3 = b/4/a*(((b/a)^2+1)*log(a/b)+(b/a)^2-1);
    C8 = 1/2*(1+v+(1-v)*(b/a)^2);
    C9 = b/a*((1+v)/2*log(a/b)+(1-v)/4*(1-(b/a)^2));
    %L3 = r0/4/a*(((r0/a)^2+1)*log(a/r0)+(r0/a)^2-1); % for case 1L
    %L9 = r0/a*((1+v)/2*log(a/r0)+(1-v)/4*(1-(r0/a)^2)); % for case 1L
    L11 = 1/64*(1+4*(r0/a)^2-5*(r0/a)^4-4*(r0/a)^2*(2+(r0/a)^2)*log(a/r0)); % for end deflection only
    L17 = 1/4*(1-(1-v)/4*(1-(r0/a)^4)-(r0/a)^2*(1+(1+v)*log(a/r0)));
    F2 = 1/4 * (1 - (b./r).^2.*(1+2*log(r/b)));
    F3 = b./(4*r) .* ( ( (b./r).^2 + 1 ).*log(r/b) + (b./r).^2 - 1);
    F8 = 1/2 * (1+v+(1-v)*(b./r).^2);
    F9 = b./r .* (1/2*(1+v)*(log(r/b)) + 1/4*(1-v)*(1-(b./r).^2));

    bracket = zeros(size(r));
    bracket(r > r0) = 1;

    ratio = r0./r;
    G11 = 1/64 * (1 + 4*ratio.^2 - 5*ratio.^4 - 4*ratio.^2.*(2+ratio.^2).*log(1./ratio) ) .* bracket;
    G17 = 1/4 * (1 - ((1-v)/4)*(1-ratio.^4) - (ratio).^2.*(1+(1+v)*log(1./ratio))) .* bracket;

    Mrb = -q*a^2/C8 * (C9*(a^2-r0^2)/(2*a*b) - L17);
    Qb = q/2/b * (a^2 - r0^2);
    E = 0;% fixme 
    D = 0; % fixme E*h^3/12/(1-v^2);

    y_over_D = Mrb * r.^2 .* F2 ...
              + Qb * r.^3 .* F3 ...
              - q  * r.^4 .* G11;
    w_nondim = y_over_D * 2*pi/(P*a^2);
    tilt_angle = 0; % fixme use equation for theta on p463

    Mr = Mrb*F8 + Qb*r.*F9 - q*r.^2.*G17;
    Mt = tilt_angle*D*(1-nu^2)./r + nu*Mr;

    Mr_nondim = Mr * 2*pi/P;
    Mt_nondim = Mt * 2*pi/P;
end

function [w_nondim,Mr_nondim,abcd] = concentrated_plate_nondim(lam,nu,theta,rho,N)
    % lam: aspect ratio: b/a = inner radius/outer radius
    % nu: poisson ratio
    % theta: angular coordinate
    % rho: nondim radial coordinate
    % N: number of terms to compute in infinite sum

    assert(isscalar(lam))
    assert(isscalar(nu))
    assert(isscalar(N))
    % vector inputs only allowed for theta and rho

    [RHO,n] = meshgrid(rho,2:N);
    d0 = -1/4;
    c0_num = (1 + 2*log(lam)) * (1+nu) - (3+nu);
    c0_den = lam^(-1)*(1+nu) + lam*(1-nu);
    c0 = lam/4 * c0_num / c0_den;
    b0 = c0/2 * (1-nu)/(1+nu) + (3+nu)/(8*(1+nu));
    a0 = -lam^2*b0 - c0*log(lam) - d0*lam^2*log(lam);
    
    d1 = 1/2;
    b1 = -1/4 * (1+nu+lam^2*(1-nu)) / (3+nu+lam^4*(1-nu));
    c1 = -b1*(3+nu)/(1-nu) - 1/4 * (1+nu)/(1-nu);
    a1 = -b1*lam^2 - c1*lam^(-2) - d1*log(lam);
    
    A = (3+nu)/(1-nu);
    B = (1-lam^2)^2 * (n.^2-1) + (lam.^(-2*n+2)+A) .* (lam.^(2*n+2)+A);
    
    dn_num = (1-lam^2)*(n-1) + lam.^( 2*n+2) + A;
    bn_num = (1-lam^2)*(n+1) - lam.^(-2*n+2) - A;
    dn_denom = B.*n.*(n-1)*(1-nu);
    bn_denom = B.*n.*(n+1)*(1-nu);
    dn = -dn_num ./ dn_denom;
    bn =  bn_num ./ bn_denom;
    an = -lam^2 * (bn .* (n+1)./n + dn .* lam.^(-2*n)./n);
    cn = -lam^2 * (dn .* (n-1)./n - bn .* lam.^( 2*n)./n);
    
    if length(rho)==1
    abcd  = [a0 b0 c0 d0;
             a1 b1 c1 d1;
             an bn cn dn];
    end
    
    Rho_zero_terms = a0 + b0 * rho.^2 + c0 * log(rho) + d0 * rho.^2 .* log(rho);
    Rho_one_terms = a1*rho + b1 * rho.^3 + c1 * rho.^(-1) + d1 * rho .* log(rho);

    
    Rho_two_to_N_terms = an .* RHO.^n + bn .* RHO.^(n+2) + cn .* RHO.^(-n) + dn .* RHO.^(-n+2);
    
    w_nondim = Rho_zero_terms' + Rho_one_terms' * cos(theta) + Rho_two_to_N_terms' * cos(n(:,1)*theta);
    
    e0_d0_coeff = 3 + nu + 2*(1+nu)*log(rho);
    e0 = 2*b0*(1+nu) - c0*rho.^-2*(1-nu) + d0*(e0_d0_coeff);
    e1 = 2*b1*rho*(3+nu) + 2*c1*rho.^-3*(1-nu) + d1*rho.^-1*(1+nu);
    
    en_a_term = an .* RHO.^(n-2) .* n .* (n-1) * (1 - nu);
    en_b_term = bn .* RHO.^n .* (n+1) .* (n + 2 - nu*(n-2));
    en_c_term = cn .* RHO.^(-n-2) .* n .* (n+1) * (1-nu);
    en_d_term = dn .* RHO.^-n .* (n-1) .* (n-2 - nu*(n+2));
    en = en_a_term + en_b_term + en_c_term + en_d_term;
    Mr_nondim = -1 * (e0' + e1'*cos(theta) + en'*cos(n(:,1)*theta));
end

function sigma_vm = damping_plate_structures_old(F_heave, D_d, D_s,P_hydrostatic,t_d,A_dt,...
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

