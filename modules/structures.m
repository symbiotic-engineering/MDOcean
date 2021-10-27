
function [B,FOS,GM] = structures(V_d, m_tot, F_hydro_heave, F_surge, F_ptrain, ...
                        M, h, rho_w, g, sigma_y, A_c, A_lat_sub, r_over_t, I, E, t_r, t_f)

%% Buoyancy Calculations
Fb = rho_w * V_d * g;
Fg = m_tot * g;
B = Fb / Fg;

%% Metacentric Height Calculatons
KB = (t_f/2 + h + t_r)/2;	% center of buoyancy above the keel
KG =  t_f/2 + h + t_r;  	% center of gravity above the keel
I_m = pi/64 * D_sft^4;      % second moment of area of the water plane area
BM = I_m / V_d;             % V_d is the submerged/displaced volume
GM = KB + BM - KG;          % Metacentric Height

%% Stress calculations
depth = h; % estimate for now

F_axial = max([F_hydro_heave, F_ptrain]);
P_hydrostatic = rho_w * g * depth;
sigma_surge = F_surge ./ A_lat_sub;

sigma_rr = P_hydrostatic + abs(sigma_surge);% radial compression
sigma_tt = P_hydrostatic * r_over_t;        % hoop stress
sigma_zz = F_axial ./ A_c;                  % axial compression
sigma_rt = sigma_surge;                     % shear
sigma_tz = [0 0 0];
sigma_zr = [0 0 0];

sigma = zeros(3,3,3);
for i=1:3
sigma(:,:,i) = [sigma_rr(i) sigma_rt(i) sigma_zr(i);
                sigma_rt(i) sigma_tt(i) sigma_tz(i);
                sigma_zr(i) sigma_tz(i) sigma_zz(i)];
end

% assume ductile material for now - need to use mohr's circle for concrete
sigma_vm = von_mises(sigma_rr, sigma_tt, sigma_zz, sigma_rt, sigma_tz, sigma_zr);

%% Buckling calculation
K = .5; % fixed-fixed - is this correct?
L = h;
F_buckling = pi^2 * E(M) * I / (K*L)^2;

%% Factor of Safety (FOS) Calculations
FOS_yield = sigma_y(M) ./ sigma_vm;
FOS_buckling = F_buckling / F_axial;
FOS = min([FOS_yield, FOS_buckling]);

end

function s_vm = von_mises(s_11, s_22, s_33, s_12, s_23, s_31)

principal_term = 1/2 * ( (s_11 - s_22).^2 + (s_22 - s_33).^2 + (s_33 - s_11).^2 );
shear_term = 3 * (s_12.^2 + s_23.^2 + s_31.^2);

s_vm = sqrt( principal_term + shear_term );

end

