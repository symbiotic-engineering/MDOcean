N = 50;
syms r theta z real
syms R_1n_1(n) R_1n_2(n) R_2n_2(n) Z_n_i1(n) Z_n_i2(n) Lambda_k(k) N_k(k) Z_k_e(k)
syms h g A omega m_k a1 a2 d1 d2 m0 real positive
syms n k real
syms C_1n_1(n) C_1n_2(n) C_2n_1(n) C_2n_2(n) B_k(k) %[N+1 1]

% eq 4
lambda_n1 = n*pi/(h-d1);
lambda_n2 = n*pi/(h-d2);

% eq 5
phi_p_i1 = 1/(2*(h-d1)) * ((z+h)^2 - r^2/2);
phi_p_i2 = 1/(2*(h-d2)) * ((z+h)^2 - r^2/2);

% eq 7
R_1n_1(n) = piecewise(n==0, 1/2, n>=1, besseli(0,lambda_n1*r)/besseli(0,lambda_n1*a2));
R_1n_2(n) = piecewise(n==0, 1/2, n>=1, besseli(0,lambda_n2*r)/besseli(0,lambda_n2*a2));

% eq 8
R_2n_1 = sym(0);
R_2n_2(n) = piecewise(n==0, 1/2*log(r/a2), n>=1, besselk(0,lambda_n2*r)/besselk(0,lambda_n2*a2));

% eq 9
Z_n_i1(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n1*(z+h)));
Z_n_i2(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n2*(z+h)));

% eq 6
% R_1n_1 = subs(R_1n_1,n,(0:N).');
% R_1n_2 = subs(R_1n_2,n,(0:N).');
% R_2n_2 = subs(R_2n_2,n,(0:N).');
% Z_n_i1 = subs(Z_n_i1,n,(0:N).');
% Z_n_i2 = subs(Z_n_i2,n,(0:N).');

sum_arg_1 = (C_1n_1 * R_1n_1 + C_2n_1 * R_2n_1) * Z_n_i1;
sum_arg_2 = (C_1n_2 * R_1n_2 + C_2n_2 * R_2n_2) * Z_n_i2;
phi_h_i1 = symsum(sum_arg_1, n, 0, N);
phi_h_i2 = symsum(sum_arg_2, n, 0, N);
% phi_h_i1 = sum(sum_arg_1);
% phi_h_i2 = sum(sum_arg_2);

phi_1 = phi_h_i1 + phi_p_i1;
phi_2 = phi_h_i2 + phi_p_i2;

% pretty(phi_1)
% pretty(phi_2)

% eq 13
Lambda_k(k) = piecewise(k==0, besselh(0,1,m0*r)/besselh(0,1,m0*a2), k>=1, besselk(0,m_k*r)/besselk(0,m_k*a2));

% eq 2.34 in analytical methods book
N_k(k) = piecewise(n==0, 1/2*(1+sinh(2*m0*h)/(2*m0*h)), n==1, 1/2*(1+sinh(2*m_k*h)/(2*m_k*h)));

% eq 14
Z_k_e(k) = piecewise(k==0, 1/sqrt(N_k) * cosh(m0 * (z+h)), k>=1, 1/sqrt(N_k) * cosh(m_k * (z+h)));

% eq 12
% sym_arg_e = B_k .* subs(Lambda_k * Z_k_e,k,(0:N).');
sym_arg_e = B_k * Lambda_k * Z_k_e;
phi_e = symsum(sym_arg_e,k,0,N);
%phi_e = sum(sym_arg_e);

% potential matching (total)
phi_1_a1 = subs(phi_1, r, a1);
phi_2_a1 = subs(phi_2, r, a1);
match_12_potential = phi_1_a1 == phi_2_a1;

phi_2_a2 = subs(phi_2, r, a2);
phi_e_a2 = subs(phi_e, r, a2);
match_2e_potential = phi_2_a2 == phi_e_a2;

% velocity matching (total)
dphi_1_dr = diff(phi_1, r);
dphi_2_dr = diff(phi_2, r);
dphi_e_dr = diff(phi_e, r);

dphi_1_dr_a1 = subs(dphi_1_dr,r,a1);
dphi_2_dr_a1 = subs(dphi_2_dr,r,a1);
match_12_velocity = dphi_1_dr_a1 == dphi_2_dr_a1;

dphi_2_dr_a2 = subs(dphi_2_dr, r, a2);
dphi_e_dr_a2 = subs(dphi_e_dr, r, a2);
match_2e_velocity = dphi_2_dr_a2 == dphi_e_dr_a2;

% body boundary condition (particular only)
dphi_p_1_dz = diff(phi_p_i1, z);
dphi_p_2_dz = diff(phi_p_i2, z);
dphi_p_1_dz_d1 = subs(dphi_p_1_dz, z, -d1);
dphi_p_2_dz_d2 = subs(dphi_p_2_dz, z, -d2);

body_1_velocity = dphi_p_1_dz_d1 == 1;
body_2_velocity = dphi_p_2_dz_d2 == 1;

eqns = [match_12_potential, match_2e_potential ...
        match_12_velocity,  match_2e_velocity ...
        body_1_velocity,    body_2_velocity];

B_n(n) = subs(B_k, k, n);
unknowns = [C_1n_1(0:N) C_1n_2(0:N) C_2n_1(0:N) C_2n_2(0:N) B_n(0:N)];

syms C_1n_1_const C_1n_2_const C_2n_1_const C_2n_2_const B_k_const [N+1 1]
unknowns_const = [C_1n_1_const; C_1n_2_const; C_2n_1_const; C_2n_2_const; B_k_const];
eqns = subs(eqns, unknowns, unknowns_const');

solve(eqns, unknowns_const)
