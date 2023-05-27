function [phi, force,added_mass,damping] = get_phi_force(eqns, unknowns_const, N, ...
    R_1n_1, R_1n_2, R_2n_1, R_2n_2, Z_n_i1, Z_n_i2, Lambda_k, Z_k_e, ...
    phi_p_i1, phi_p_i2, r, z, spatial_res, plot_phi, ...
    h_num, m_k_num, a1_num, a2_num, d1_num, d2_num, m0_num)

syms m_k_const [1 N]
eqns = subs(eqns,m_k_const,m_k_num);

syms h a1 a2 d1 d2 m0 real positive
params = {h a1 a2 d1 d2 m0};
params_num = {h_num a1_num a2_num d1_num d2_num m0_num};
eqns = subs(eqns,params,params_num);

solns = vpasolve(eqns, unknowns_const);

fns = fieldnames(solns);

C_1n_1s = ones(N+1,1);
C_2n_1s = zeros(N+1,1); % because R_2n_1s = 0
C_1n_2s = ones(N+1,1);
C_2n_2s = ones(N+1,1);
B_ks = ones(N+1,1);
index = 1;

for k = 1:length(fns)
    if index > N+1
        index = 1;
    end
    if isempty(solns.(fns{k}))
        error('Linear solver could not find a solution')
    else
        if (k<=N+1)
            C_1n_1s(index) = solns.(fns{k});
        end
        if (N+1<k) && (k<=2*(N + 1))
            C_1n_2s(index) = solns.(fns{k});
        end
        if (k>2*(N+1)) && (k<=3*(N+1))
            C_2n_2s(index) = solns.(fns{k});
        end
        if (k>3*(N+1))
            B_ks(index) = solns.(fns{k});
        end
    end
    index = index + 1;
end

phi_h_n_i1_solns_all = (C_1n_1s' .* R_1n_1(0:N) + C_2n_1s' .* R_2n_1(0:N)) .* Z_n_i1(0:N);
phi_h_n_i2_solns_all = (C_1n_2s' .* R_1n_2(0:N) + C_2n_2s' .* R_2n_2(0:N)) .* Z_n_i2(0:N);

phi_e_k = B_ks' .* Lambda_k .* Z_k_e;

% summing all to get the potential
phi_h_i1 = sum(phi_h_n_i1_solns_all,2);
phi_h_i2 = sum(phi_h_n_i2_solns_all,2);

phi_e_sym = sum(phi_e_k);

% sub in geometry
phi_h1 = subs(phi_h_i1,params,params_num);
phi_h2 = subs(phi_h_i2,params,params_num);
phi_p1 = subs(phi_p_i1,params,params_num);
phi_p2 = subs(phi_p_i2,params,params_num);
phi_e = subs(phi_e_sym,params,params_num);
phi_e = subs(phi_e,m_k_const,m_k_num);

% 36 and 37
phi_1 = phi_h1 + phi_p1;
phi_2 = phi_h2 + phi_p2;

% velocity
v_1r = diff(phi_1,r);
v_1z = diff(phi_1,z);
v_2r = diff(phi_2,r);
v_2z = diff(phi_2,z);
v_er = diff(phi_e,r);
v_ez = diff(phi_e,z);

% sub in spatial coordinates
r_vec = linspace(2*a2_num/spatial_res,2*a2_num,spatial_res);
r_vec = sort([r_vec a1_num a2_num]);
z_vec = linspace(-h_num,0,spatial_res);
[R,Z] = meshgrid(r_vec,z_vec);

phi1 = double(subs(phi_1,{r,z},{R,Z}));
phi2 = double(subs(phi_2,{r,z},{R,Z}));
phi1h = double(subs(phi_h1,{r,z},{R,Z}));
phi1p = double(subs(phi_p1,{r,z},{R,Z}));
phi2h = double(subs(phi_h2,{r,z},{R,Z}));
phi2p = double(subs(phi_p2,{r,z},{R,Z}));
phie = double(subs(phi_e,{r,z},{R,Z}));

v1r = double(subs(v_1r,{r,z},{R,Z}));
v1z = double(subs(v_1z,{r,z},{R,Z}));
v2r = double(subs(v_2r,{r,z},{R,Z}));
v2z = double(subs(v_2z,{r,z},{R,Z}));
ver = double(subs(v_er,{r,z},{R,Z}));
vez = double(subs(v_ez,{r,z},{R,Z}));

%calculate added mass and damping
dphi_i_dz = v1z + v2z+ vez; % matrix
disp(dphi_i_dz)
phi_i = phi1 + phi2 ;%+ phie;
rho = 1023;
%region by region
added_mass1 = real(1023*(h^3)*2*pi*(int(r,0,a1)+int(r,0,a2)))*phi_i*dphi_i_dz); %+ int(r,0,a2)) %dimension error

damping =   imag(1023*(h^3)*2*pi*int(r,0,a1)*phi_i*dphi_i_dz);

%non dimensionalizing it
added_mass = added_mass / (rho*pi*a1)

disp(added_mass)



% assemble total phi based on phi in each region
regione = R > a2_num;
region1 = R <= a1_num & Z < -d1_num;
region2 = R > a1_num & R <= a2_num & Z < -d2_num;

phi = NaN(size(R));
phi(region1) = phi1(region1);
phi(region2) = phi2(region2);
phi(regione) = phie(regione);

phiH = NaN(size(R));
phiH(region1) = phi1h(region1);
phiH(region2) = phi2h(region2);

phiP = NaN(size(R));
phiP(region1) = phi1p(region1);
phiP(region2) = phi2p(region2);

v_r = NaN(size(R));
v_r(region1) = v1r(region1);
v_r(region2) = v2r(region2);
v_r(regione) = ver(regione);

v_z = NaN(size(R));
v_z(region1) = v1z(region1);
v_z(region2) = v2z(region2);
v_z(regione) = vez(regione);

if plot_phi
    region_body = ~region1 & ~region2 & ~regione;
    plot_potential(phi,R,Z,region_body,'Total Potential');
    plot_potential(phiH,R,Z,region_body,'Homogeneous Potential');
    plot_potential(phiP,R,Z,region_body,'Particular Potential');
    plot_potential(v_r,R,Z,region_body,'Radial Velocity')
    plot_potential(v_z,R,Z,region_body,'Vertical Velocity')
    plot_velocity(real(v_r),real(v_z),R,Z);

    plot_matching(phi1,phi2,phie,a1_num,a2_num,R,Z,'\phi')
    plot_matching(v1r,v2r,ver,a1_num,a2_num,R,Z,'v_r')
end

% force
force_2_over_rho_w_eiwt = 1i * int( int(r*(phi_h_i2+phi_p_i2),a1,a2), 0, 2*pi );

%symvar(force_2_over_rho_w_eiwt)

force = real(double(vpa(subs(force_2_over_rho_w_eiwt,{h,a1,a2,d2},{h_num,a1_num,a2_num,d2_num}))));
end