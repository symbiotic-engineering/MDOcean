function [V_d, V_m, m_tot, m_float, h] = geometry(D_i, D_sft, t_sft, L_sf, ...
                                t_sf, t_sfb, t_vc, D_or, t_r, rho_m, M)

%D_sft-diameter of the top surface float plate
%t_sft-thickness of the top surface float plate
%L_sf- length of the surface float column (rectangle)
%W_sf-width of the surface float column (rectangle)
%t_sf-thickness of the surface float column
%t_sfb-thickness of the bottom surface float plate
                            
%% Surface float
%Volume of the surface float: (for displacement purposes)
t_f = D_sft / 4; % scaling law to keep same proportions as RM3
V_sf_d = 1/2 * pi * ((D_sft/2)^2 - (D_i/2)^2) * t_f; % 1/2 factor to account for waterline (not fully submerged)
%Volume of the surface float: (for material purposes)
W_sf = pi * D_i; % circumference
V_sf_m = (pi * D_sft/2 * t_sft) + L_sf * W_sf * t_sf + (pi * D_sft/2 * t_sfb);
m_float = V_sf_m * rho_m(M);

%% Vertical column
h = D_i / 7; % scaling law to keep same proportions as RM3
%Volume of the vertical column: (for displacement purposes)
V_vc_d = pi * (D_i/2)^2 * h;
%Volume of the vertical column: (for material purposes)
D_ivc = D_i - t_vc;
V_vc_m = pi * ((D_i/2)^2-(D_ivc/2)^2) * h;

%% Reaction plate
%Volume of the reaction plate (for displacement and material purposes):
V_rp = pi * ((D_or/2)^2 - (D_i/2)^2) * t_r;

%% Totals
V_d = V_sf_d + V_vc_d + V_rp; % Total Volume displaced
V_m = V_sf_m + V_vc_m + V_rp; % Total Volume of material
m_tot = V_m * rho_m(M);       % Total mass of material

end

