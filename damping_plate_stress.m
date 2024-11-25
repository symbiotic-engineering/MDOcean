    syms D_d A_dt theta_dt L_dt h_d D_s t_d
    syms A_c [1,3]
    syms P_hydrostatic [1,3]
    syms E [1,3]

    %What codey thing dictates which material we're using so we don't have
    %to hardcode this?
    E = E(1);

    %calculate the deflection of the damping plate
    r = D_d/2;
    r_bending = r-(D_s/2);
    F_water = P_hydrostatic(2) * A_c(3)/2;
    syms F_supportsym xsym
    Msym = -(F_water+F_supportsym*sin(theta_dt))*xsym + (F_water*xsym^2/(2*r_bending)) + ...
        (F_water*r_bending/2) + (F_supportsym*r_bending*sin(theta_dt));
    Isym = (1/6) * sqrt(r_bending^2 - xsym^2) * (h_d^3 - t_d^3);
    y_double_prime = Msym/(E*Isym);
    y_prime = int(y_double_prime,xsym);
    y = int(y_prime,xsym);
    y_tip = subs(y,xsym,r_bending-0.5);
    eqn = F_supportsym == y_tip * E * A_dt / L_dt;
    F_support = solve(eqn,F_supportsym);

    I_d = (1/12) * r_bending * (h_d^3 - t_d^3);
    %F_support = 3*F_water*r_bending^3*A_dt/(8*sin(theta_dt)*(3*L_dt*I_d - (r_bending^3*A_dt)));

    % x1 = linspace(0,r-t_d,1500);
    % A = 2*sqrt(r^2-x1.^2) * h_d - (2*sqrt(r^2-x1.^2) * t_d);
    % sigma_axial = F_support * cos(theta_dt) ./ A;
    % v = abs(-F_water - (F_support*(sin(theta_dt))) + (F_water*x1/r));
    % shear = max(v./A);

    x = linspace(0,r_bending-0.01,1500);
    x1 = linspace(0,r-t_d,1500);

    A = 2*sqrt(r_bending^2-x.^2) * h_d - (2*sqrt(r_bending^2-x.^2) * t_d);
    sigma_axial = F_support * cos(theta_dt) ./ A;
    v = -F_water - (F_support*(sin(theta_dt))) + (F_water*x1/r_bending);
    shear = max(v./A);

    M_x = -(F_water+F_support*sin(theta_dt))*x + (F_water.*x.^2/(2*r_bending)) + ...
        (F_water*r_bending/2) + (F_support*r_bending*sin(theta_dt));
    I_x = (1/12) * 2*sqrt(r_bending^2-x.^2) * h_d^3 - (((1/12) * 2*sqrt(r_bending^2-x.^2) * t_d^3));
    sigma_bending = abs(M_x .* h_d./(2.*I_x));
    sigma_xx = max(sigma_bending+sigma_axial);

    %debugging
    % figure
    % hold on
    % title("Bending and Axial Stress vs. Radial Position")
    % plot(x,sigma_bending)
    % plot(x,sigma_axial)
    % legend("bending","axial")
    % xlabel("x")
    % ylabel("stress")
    % hold off
    % 
    % figure
    % hold on
    % title("Shear and Bending Moment Diagram")
    % plot(x,v)
    % plot(x,M_x)
    % legend("shear","bending moment")
    % xlabel("x")
    % hold off
    % 
    % figure
    % hold on
    % title("Moment of Inertia vs. Radial Position")
    % plot(x,I_x)
    % legend("moment of inertia")
    % xlabel("x")
    % ylabel("I")
    % hold off
    % 
    % figure
    % hold on
    % plot(x,sigma_axial)
    % title("Axial Stress vs. Radial Position")
    % xlabel("x")
    % ylabel("Axial Stress")
    % hold off

    ht = matlabFunction(sigma_xx, shear, 'vars', [E, D_d, A_dt, theta_dt, L_dt, h_d, D_s, t_d, A_c, P_hydrostatic])