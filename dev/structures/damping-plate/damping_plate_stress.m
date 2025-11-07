syms D_d A_dt theta_dt L_dt h_d D_s t_d F_heave
syms A_c [1,3]
syms P_hydrostatic [1,3]
syms E [1,3]
syms x x1 c1 c2

%What codey thing dictates which material we're using so we don't have
%to hardcode this?
E = E(1);

%set up the dimensions and pressures
r = D_d/2;
r_bending = r-(D_s/2);
%F_water = P_hydrostatic(2) * A_c(3)/2;
P_hydrodynamic = F_heave / A_c(3);
F_water = P_hydrodynamic * A_c(3) / 2;

%solve y" = M/EI to find the deflection of the damping plate
syms F_supportsym xsym
Msym = -(F_water+F_supportsym*sin(theta_dt))*xsym + (F_water*xsym^2/(2*r_bending)) + ...
    (F_water*r_bending/2) + (F_supportsym*r_bending*sin(theta_dt));
Isym = (1/6) * sqrt(r_bending^2 - xsym^2) * (h_d^3 - t_d^3);
y_double_prime = Msym/(E*Isym);
y_prime = int(y_double_prime,xsym) + c1;
find_c1 = 0 == subs(y_prime,xsym,0);
c1_solved = solve(find_c1,c1);
y_prime = subs(y_prime, c1, c1_solved);
y = int(y_prime,xsym) + c2;
find_c2 = 0 == subs(y, xsym, 0);
c2_solved = solve(find_c2,c2);
y = subs(y,c2,c2_solved);

%use the tip deflection (which is equal to the deflection of the tubular
%support) to find the force of the tubular support on the damping plate
y_tip = subs(y,xsym,r_bending-1);
eqn = F_supportsym == y_tip * E * A_dt / L_dt;
F_support = solve(eqn,F_supportsym);

%solve for deflection as a function of x for debugging
y = subs(y, F_supportsym, F_support);
y = subs(y,xsym,x);

%solve for the shear and axial stresses
A = 2*sqrt(r_bending^2-x^2) * h_d - (2*sqrt(r_bending^2-x^2) * t_d);
sigma_axial = F_support * cos(theta_dt) ./ A;
v = -F_water - (F_support*(sin(theta_dt))) + (F_water*x1/r_bending);
shear = max(v./A);

%solve for the bending stress
M_x = -(F_water+F_support*sin(theta_dt))*x + (F_water.*x^2/(2*r_bending)) + ...
    (F_water*r_bending/2) + (F_support*r_bending*sin(theta_dt));
I_x = (1/12) * 2*sqrt(r_bending^2-x^2) * h_d^3 - (((1/12) * 2*sqrt(r_bending^2-x^2) * t_d^3));
sigma_bending = M_x .* h_d./(2.*I_x);
sigma_xx = sigma_bending+sigma_axial;

matlabFunction(sigma_xx, shear, sigma_bending, sigma_axial, M_x,I_x, F_water, F_support, y, 'File','damping_plate_func','Vars', [E, D_d, A_dt, theta_dt, L_dt, h_d, D_s, t_d, A_c, P_hydrostatic,x,x1,F_heave]);
