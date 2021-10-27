function [F_heave,F_surge, m, D_env,V_s] = hydro(t, s, m_float, D_sft, D_i,D_or,t_r,d_WEC, N_WEC, d_farm, d_shore, rho_w, g, Hs, T, t_f)

w = 2*pi/T;        % frequency
k = w^2 / g;        % wave number
V_g = g /(2*w);     % group velocity

%Surface Float Surge Geometry Calculations
r_sf = D_sft / 2;      % surface float radius
A_wsf = pi * r_sf^2;     % surface float waterplane area
draft_sf = t_f / 2;    % surface float estimated depth
vol_sf = A_wsf * draft_sf;  % surface float submerged volume

%Vertical Column Surge Geometry Calculations
r_vc=D_i/2;% vertical column radius
A_wvc=pi*r_vc^2;%vertical column waterplane area
draft_vc=D_i*7;%vertical column depth
vol_vc=A_wvc*draft_vc;%vertical column submerged volume

%Reaction Plate Geometry Calculations
r_rp=D_or/2;%reaction plate radius
A_wrp=pi*r_rp^2;%reation plate waterplane area
draft_rp=t_r;%reaction plate depth
vol_rp=A_wrp*draft_rp;%vertical column submerged volume

%Submerged RM3 Volume
V_s=[vol_sf vol_vc vol_rp];

% hydrodynamic coefficients
A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63 * ones(1,size(s,2)); % added mass
gamma   = rho_w * g * A_wsf * ones(1,size(s,2)); % Froude Krylov / diffraction
B       = k / (4 * rho_w * g * V_g) * gamma.^2; % radiation damping

wave = Hs * sin(w * t); % in reality this is a function of d_shore

F_hydrostatic = -rho_w * g * A_wsf * s(1,:);
F_radiation = -B .* s(2,:);
F_excitation = gamma .* wave;
F_heave = F_hydrostatic + F_radiation + F_excitation;

%Surge Force Per Geometry Calculations
F_surge_sf= 2 * (Hs/2) * rho_w * g * vol_sf* (1 - exp(-k*(draft_sf)));%surface float surge force
F_surge_vc= 2 * (Hs/2) * rho_w * g * vol_vc * (1 - exp(-k*(draft_vc)));%vertical column surge force
F_surge_rp= 2 * (Hs/2) * rho_w * g * vol_rp * (1 - exp(-k*(draft_rp)));%reaction plate surge force
F_surge=[F_surge_sf F_surge_vc F_surge_rp];%Surge force Matrix [surface float, vertical column, reaction plate]
m = m_float + A;

D_env = 0.5 * ones(1,size(s,2)); % in reality this is a function of d_farm and forces
end

