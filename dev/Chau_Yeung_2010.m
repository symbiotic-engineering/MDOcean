N = 20;
syms r theta z real
syms R_1n_1(n) R_1n_2(n) R_2n_2(n) Z_n_i1(n) Z_n_i2(n) Lambda_k(k) N_k(k) Z_k_e(k)
syms h m_k a1 a2 d1 d2 m0 real positive
syms n k real positive
syms C_1n_1(n) C_1n_2(n) C_2n_1(n) C_2n_2(n) B_k(k) %[N+1 1]

% eq 4
lambda_n1(n) = n*pi/(h-d1);
lambda_n2(n) = n*pi/(h-d2);

% eq 5
phi_p_i1 = 1/(2*(h-d1)) * ((z+h)^2 - r^2/2);
phi_p_i2 = 1/(2*(h-d2)) * ((z+h)^2 - r^2/2);

% eq 7
R_1n_1(n) = piecewise(n==0, 1/2, n>=1, besseli(0,lambda_n1(n)*r)/besseli(0,lambda_n1(n)*a2));
R_1n_2(n) = piecewise(n==0, 1/2, n>=1, besseli(0,lambda_n2(n)*r)/besseli(0,lambda_n2(n)*a2));

% eq 8
R_2n_1(n) = sym(0);
R_2n_2(n) = piecewise(n==0, 1/2*log(r/a2), n>=1, besselk(0,lambda_n2(n)*r)/besselk(0,lambda_n2(n)*a2));

% eq 9
Z_n_i1(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n1(n)*(z+h)));
Z_n_i2(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n2(n)*(z+h)));

% eq 6
phi_h_n_i1(n) = (C_1n_1(n) * R_1n_1(n) + C_2n_1(n) * R_2n_1(n)) * Z_n_i1(n);
phi_h_n_i2(n) = (C_1n_2(n) * R_1n_2(n)+ C_2n_2(n) * R_2n_2(n)) * Z_n_i2(n);
%phi_h_i1 = symsum(phi_h_n_i1, n, 0, N);
%phi_h_i2 = symsum(phi_h_n_i2, n, 0, N);
% phi_h_i1 = sum(phi_h_n_i1);
% phi_h_i2 = sum(phi_h_n_i2);

%phi_1 = phi_h_i1 + phi_p_i1;
%phi_2 = phi_h_i2 + phi_p_i2;

% pretty(phi_1)
% pretty(phi_2)

% eq 13
Lambda_k(k) = piecewise(k==0, besselh(0,1,m0*r)/besselh(0,1,m0*a2), k>=1, besselk(0,m_k*r)/besselk(0,m_k*a2));

% eq 2.34 in analytical methods book
N_k(k) = piecewise(k==0, 1/2*(1+sinh(2*m0*h)/(2*m0*h)), k>=1, 1/2*(1+sinh(2*m_k*h))/(2*m_k*h));

% eq 14
Z_k_e(k) = piecewise(k==0, 1/sqrt(N_k(k)) * cosh(m0 * (z+h)), k>=1, 1/sqrt(N_k(k)) * cosh(m_k * (z+h)));

% eq 12
% sym_arg_e = B_k .* subs(Lambda_k * Z_k_e,k,(0:N).');
%phi_e_k = B_k * Lambda_k * Z_k_e;
%phi_e = symsum(phi_e_k,k,0,N);
%phi_e = sum(phi_e_k);

% potential matching (total)
% phi_1_a1 = subs(phi_1, r, a1);
% phi_2_a1 = subs(phi_2, r, a1);
%match_12_potential = phi_1_a1 == phi_2_a1;

B_n(n) = subs(B_k, k, n);
Lambda_n(n) = subs(Lambda_k, k, n);
Z_n_e(n) = subs(Z_k_e, k, n);

% equation 22 in old 1981 paper, applied to boundary 2-e
dz_2 = 1 - d2/h;
match_2e_potential = C_1n_2(n) * subs(R_1n_2(n),r,a2) + C_2n_2(n) * subs(R_2n_2(n),r,a2) == ...
    B_n(n) * subs(Lambda_n(n),r,a2) * dz_2 - int(subs(phi_p_i2,r,a2) * Z_n_i2(n), z, 0, dz_2);

% equation 22 in old 1981 paper, applied to boundary 1-2
dz_1 = 1 - d1/h;
match_12_potential = C_1n_1(n) * subs(R_1n_1(n),r,a1) == ...
    ( C_1n_2(n) * subs(R_1n_2(n),r,a1) + C_2n_2(n) * subs(R_2n_2(n),r,a1) ) * dz_1 + ...
    int(subs(phi_p_i2 - phi_p_i1,r,a1) * Z_n_i1(n), z, 0, dz_1);

%phi_2_a2 = subs(phi_2, r, a2);
%phi_e_a2 = subs(phi_e, r, a2);
%match_2e_potential = phi_2_a2 == phi_e_a2;

% velocity matching (total)

% equation 23 in old 1981 paper, applied to boundary 2-e - giving NAN for
% N=2
match_2e_velocity = B_n(n) * subs(diff(Lambda_n(n), r), r, a2) == ...
    (C_1n_2(n) * subs(diff(R_1n_2(n), r), r, a2) + C_2n_2(n) * subs(diff(R_2n_2(n), r), r, a2) ) * dz_2 + ...
    int(subs(diff(phi_p_i2,r),r,a2) * Z_n_e(n), z, 0, dz_2 );

% equation 23 in old 1981 paper, applied to boundary 1-2
match_12_velocity = C_1n_2(n) * subs(R_1n_2(n), r, a1) + C_2n_2(n) * subs(R_2n_2(n), r, a1) == ...
    C_1n_1(n) * subs( diff(R_1n_1(n),r), r,a1) * dz_1 + int( subs(diff(phi_p_i1 - phi_p_i2,r),r,a1) * Z_n_i2(n), z, 0, dz_1 );

% dphi_1_dr = diff(phi_1, r);
% dphi_2_dr = diff(phi_2, r);
% dphi_e_dr = diff(phi_e, r);
% 
% dphi_1_dr_a1 = subs(dphi_1_dr,r,a1);
% dphi_2_dr_a1 = subs(dphi_2_dr,r,a1);
% match_12_velocity = dphi_1_dr_a1 == dphi_2_dr_a1;
% 
% dphi_2_dr_a2 = subs(dphi_2_dr, r, a2);
% dphi_e_dr_a2 = subs(dphi_e_dr, r, a2);
% match_2e_velocity = dphi_2_dr_a2 == dphi_e_dr_a2;

eqns = [subs(match_12_potential,n,0:N), subs(match_2e_potential,n,0:N) ...
        subs(match_12_velocity,n,0:N),  subs(match_2e_velocity,n,0:N)];

%eqns = subs(eqns, n, 0:N);
unknowns = [C_1n_1(0:N) C_1n_2(0:N) C_2n_1(0:N) C_2n_2(0:N) B_n(0:N)];

syms C_1n_1_const C_1n_2_const C_2n_1_const C_2n_2_const B_k_const [N+1 1] real
unknowns_const = [C_1n_1_const; C_1n_2_const; C_2n_1_const; C_2n_2_const; B_k_const];
eqns = subs(eqns, unknowns, unknowns_const');

h_num = 1;
m_k_num = 1;
a1_num = .5;
a2_num = 1;
d1_num = .5;
d2_num = .25;
m0_num = 1;
params = {h m_k a1 a2 d1 d2 m0};
params_num = {h_num m_k_num a1_num a2_num d1_num d2_num m0_num};
eqns = subs(eqns,params,params_num);

solns = vpasolve(eqns, unknowns_const);
solns

% build this...
%C_1n_1_const_list = solns[0:N] 
%eq6
%phi_h_n_i1_solns = (solns.C_1n_1_const1 * R_1n_1(0) + solns.C_2n_1_const1 * R_2n_1(0)) * Z_n_i1(0);
%phi_h_n_i2_solns = (solns.C_1n_2_const2 * R_1n_2(0) + solns.C_2n_2_const2 * R_2n_2(0)) * Z_n_i2(0);

fns = fieldnames(solns);
length(fns)

C_1n_1s = ones(N+1,1);
C_2n_1s = ones(N+1,1);
C_1n_2s = ones(N+1,1);
C_2n_2s = ones(N+1,1);
B_ks = ones(N+1,1);
index = 1;

for k = 1:length(fns)
    if index > N+1
        index = 1;
    end
   
    if (k<=N+1)
        C_1n_1s(index) = solns.(fns{k});
    end
    if (N+1<k) && (k<=2*(N + 1))
        C_1n_2s(index) = solns.(fns{k});
       
    end
    if (k>2*(N+1)) && (k<=3*(N +1) )
         C_2n_1s(index) = solns.(fns{k});
         
    end

    if (k>3*(N+1)) && (k<=4*(N+1))
         C_2n_2s(index) = solns.(fns{k});
       
    end

    if (k>4*(N+1)) 
         B_ks(index) = solns.(fns{k});
    end
    index = index + 1;


end



phi_h_n_i1_solns_all = (C_1n_1s' .* R_1n_1(0:N) + C_2n_1s' .* R_2n_1(0:N)) .* Z_n_i1(0:N); 
phi_h_n_i2_solns_all = (C_1n_2s' .* R_1n_2(0:N) + C_2n_2s' .* R_2n_2(0:N)) .* Z_n_i2(0:N);

phi_e_k = B_ks' .* Lambda_k(0:N) .* Z_k_e(0:N);


% summing all to get the potential 
phi_h_i1 = sum(phi_h_n_i1_solns_all,2);
phi_h_i2 = sum(phi_h_n_i2_solns_all,2);

phi_e = sum(phi_e_k);

% 36 and 37
phi_1 = phi_h_i1 + phi_p_i1;
phi_2 = phi_h_i2 + phi_p_i2;

phi_1 = subs(phi_1,params,params_num);
phi_2 = subs(phi_2,params,params_num);
phi_e = subs(phi_e,params,params_num);


%pretty(phi_1)
%pretty(phi_2)

res = 30;
r_vec = linspace(0,2*a2_num,res);
z_vec = linspace(0,h_num,res);
[R,Z] = meshgrid(r_vec,z_vec);

phi1 = double(subs(phi_1,{r,z},{R,Z}));
phi2 = double(subs(phi_2,{r,z},{R,Z}));
phie = double(subs(phi_e,{r,z},{R,Z}));

regione = R > a2_num;
region1 = R < a1_num & Z < (h_num - d1_num);
region2 = R > a1_num & R <= a2_num & Z < (h_num - d2_num);

phi = NaN(size(R));
phi(region1) = phi1(region1);
phi(region2) = phi2(region2);
phi(regione) = phie(regione);
%%
figure
levels = [linspace(0,2,5) 3:9];
subplot 121
[c,h] = contourf(R,Z,real(phi),levels);
clabel(c,h)
xlabel('R')
ylabel('Z')
title('Velocity Potential - Real')
colorbar

imag_phi = imag(phi);
imag_phi(~region1 & ~region2 & ~regione) = NaN;

subplot 122
[c,h] = contourf(R,Z,imag_phi);
clabel(c,h)
xlabel('R')
ylabel('Z')
title('Velocity Potential - Imaginary')
colorbar

