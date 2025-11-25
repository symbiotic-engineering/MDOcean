function  figs = viz_damping_plate()

    %% Plot deflection and bending moment for RM3 geometry, just point loads
    
    b = 3; % inner radius
    a = 15; % outer radius
    h = 0.0254; % thickness
    h_mult = 10; % multiplier on h to get equivalent plate thickness including stiffeners
    h_eq = h * h_mult;
    
    F_heave = 8.5e6; % force distributed on bottom
    P = F_heave/4; % force per support
    
    E = 210e9;
    sigma_y = 248e6;
    nu = 0.3;
    
    D = E * h_eq^3 / (12*(1-nu^2));
    lam = b/a;
    theta = 0:pi/32:pi;
    rho = linspace(lam,1,15);
    N = 100;
    
    [w_nondim,Mr_nondim] = concentrated_plate_nondim(lam,nu,theta,rho,N);
    
    Mr_dim = Mr_nondim * P/(2*pi);
    sigma_dim = Mr_nondim * 3/pi * P/h^2;
    FOS = sigma_y ./ sigma_dim;
    FOS(abs(FOS)>3) = sign(FOS(abs(FOS)>3)) * 3;
    w_dim = w_nondim * P*a^2/(2*pi*D);
    
    w_normalized = w_nondim / max(abs(w_nondim),[],'all');
    Mr_normalized = Mr_nondim / max(abs(Mr_nondim),[],'all');
    
    idx_max = 1; % max bending moment occurs at inner edge
    figs_Mr_norm = make_four_loads_plot(theta,rho,Mr_normalized,idx_max,'Normalized bending moment (-)',true);
    figs_Mr_dim  = make_four_loads_plot(theta,rho,Mr_dim/1e6,idx_max,'Bending moment (MNm)',false);
    figs_sigma   = make_four_loads_plot(theta,rho,sigma_dim/1e6,idx_max,'Bending stress (MPa)',false);
    figs_FOS     = make_four_loads_plot(theta,rho,FOS,idx_max,'FOS (-)',true);
    
    idx_max = length(rho); % max deflection occurs at outer edge
    %figs_w_norm = make_four_loads_plot(theta,rho,w_normalized,idx_max,'Normalized deflection (-)');
    figs_w_dim   = make_four_loads_plot(theta,rho,w_dim,idx_max,'Deflection (m)',false);
    
    %% Superpose distributed loading
    [w_nondim,Mr_nondim] = calc_plate_superposed(a,b,F_heave,nu,theta,rho,N);
    w_dim = w_nondim * P*a^2/(2*pi*D);
    Mr_dim = Mr_nondim * P/(2*pi);
    sigma_dim = Mr_nondim * 3/pi * P/h^2;
    
    idx_max = 1; % max bending moment occurs at inner edge
    figs_Mr_dim_super = make_four_loads_plot(theta,rho,Mr_dim/1e6,idx_max,'Bending moment (MNm)',false);
    figs_sigma_super  = make_four_loads_plot(theta,rho,sigma_dim/1e6,idx_max,'Bending stress (MPa)',false);
    %figs_FOS_super   = make_four_loads_plot(theta,rho,FOS,idx_max,'FOS (-)',true);
    
    idx_max = length(rho); % max deflection occurs at outer edge
    %figs_w_norm_super = make_four_loads_plot(theta,rho,w_normalized,idx_max,'Normalized deflection (-)');
    figs_w_dim_super = make_four_loads_plot(theta,rho,w_dim,idx_max,'Deflection (m)',false);
    
    %% sweep aspect ratio
    f_aspect_sweep = plate_aspect_ratio_sweep(theta,a,F_heave,b,N,nu);
    
    %% compile figures
    % 5 runs concentrated + 3 runs superposed = 8 runs total. 8 * 3 figs each = 24 figs
    figs_four_loads_conc = [figs_Mr_norm, figs_Mr_dim, figs_sigma, figs_FOS, figs_w_dim];
    figs_four_loads_super = [figs_Mr_dim_super, figs_sigma_super, figs_w_dim_super];
    figs_four_loads = [figs_four_loads_conc, figs_four_loads_super];
    
    figs = [figs_four_loads,f_aspect_sweep]; % 25 figs total

end

%% helper functions
function f = plate_aspect_ratio_sweep(theta,a,F_heave,b_nom,N,nu)

    % Keep outer radius (a) fixed and sweep inner radius (b). Plot max stress
    % and deflection.
    
    % inner radius
    b = [.75:.25:2.5 3:14 14.5:.05:14.95];
    
    for i = 1 : length(b)
        rho = [b(i)/a 1];
        [w_temp,Mr_temp] = calc_plate_superposed(a,b(i),F_heave,nu,theta,rho,N);
        w_nondim(i) = abs(w_temp(end,1));
        Mr_nondim(i) = abs(Mr_temp(1,1));
    end
    
    
    f = figure;
    plot(b/a, w_nondim/w_nondim(b==b_nom))
    hold on
    plot(b/a, Mr_nondim/Mr_nondim(b==b_nom))
    plot(b_nom/a*[1 1],[0 1], 'k--',[0 b_nom/a], [1 1], 'k--')
    xlabel('Plate inner radius to outer radius ratio')
    ylabel('Normalized maximum')
    legend('Deflection','Stress','Nominal')
    grid on
    improvePlot

end

function [w_nondim,Mr_nondim] = calc_plate_superposed(a,b,F_heave,nu,theta,rho,N)
    % fixme: this assumes that F_tube = F_heave/4, ie that all force goes into
    % the tubes and none goes into the bottom of the spar, which is wrong.
    % Should actually use compatibility and tube stiffness to get F_tube, as in 
    % damping_plate_structures() in structures.m.
    %     lam = b/a;
    %     [w_nondim_conc,Mr_nondim_conc] = calc_plate(lam,nu,theta,rho,N);
    %     [w_nondim_dist,Mr_nondim_dist] = plate_distributed(a,b,F_heave,nu,rho);
    %     w_nondim = w_nondim_conc   + w_nondim_dist';
    %     Mr_nondim = Mr_nondim_conc + Mr_nondim_dist';
    
    vb = var_bounds();

    D_d = 2*a;
    D_s = 2*b;
    [P_hydrostatic,A_dt,theta_dt,h_d,A_c] = deal(0); % not used, so set to 0
    t_d = vb.t_d_nom * 1e-3;
    L_dt = (D_d - D_s) / (2*cos(theta_d_tu)); % formula copy-pasted from geometry.m
    h_stiff     = p.h_over_h1_stiff_d * vb.h1_stiff_d;
    width_stiff = p.w_over_h1_stiff_d * vb.h1_stiff_d;
    
    [sigma_vm,delta_max] = damping_plate_structures(F_heave, D_d, D_s,P_hydrostatic,t_d,A_dt,...
                                                theta_dt, L_dt, h_d, A_c, p.E(1), nu, h_stiff, width_stiff,...
                                                p.D_d_tu, p.t_d_tu, p.num_terms_plate, p.radial_mesh_plate, ...
                                                p.num_stiff_d, true);


    [Mr_over_F_heave, ...
          Mt_over_F_heave, ...
          Q_nondim, ...
          delta_over_a2Fheave_2piDeq, ...
          F_tube_over_F_heave]     = combined_plate_nondim(rho, b, a, nu, N, plate_tube_K_ratio);

end
