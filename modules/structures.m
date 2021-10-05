
function [B,FOS] = structures(V_d, m_tot, F_max, M, D_i, D_sft, h, D_or, rho_w, g, sigma_y)

%% Buoyancy Calculations
Fb = rho_w * V_d * g;
Fg = m_tot * g;
B = Fb / Fg;

%% Factor of Safety (FOS) Calculations
% Surface float
A_sf = pi*((D_sft/2)^2-(D_i/2)^2);
sigma_sf = F_max / A_sf;

%Vertical column: 
A_vc = 2*pi*(D_i/2)^2+h*(2*pi*(D_i/2));
sigma_vc = F_max / A_vc;

%Design stress of the reaction plate: 
A_rp = pi*((D_or/2)^2 - (D_i/2)^2);
sigma_rp = F_max / A_rp;%<--final equation

% Overall
max_stress = max([sigma_sf, sigma_vc, sigma_rp]);
FOS = sigma_y(M) / max_stress;

end

