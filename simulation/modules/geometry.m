function [V_d, V_m, m_tot, m_float, h, A_c, ...
         A_lat_sub, r_over_t, I, draft, B_f, B_s, GM] = geometry(D_s, D_f, T_f, ...
                                                  t_ft, t_fr, t_fb, t_sr, ...
                                                  rho_m, M, rho_w)

%% Variable Definitions
% D: diameter
% T: draft
% h: height
% t: material thickness
% _f: float
% _s: spar
% _d: damping plate

%               _____                               -
%               |   |                               |
%               |   |                               |
%            _____Df_____                       -   |
%           |           |                       |   |
% ----------|           |---------  -   -       hf  |
%           |           |           Tf  |       |   |
%           _____________           -   |       -   |
%               |   |                   |           |
%               |   |                   |           hs
%               |Ds |                   Ts          |
%               |   |                   |           |
%               |   |                   |           |
%               |   |                   |           |
%               |   |                   |           |
%       _________Dd__________           -   -       -
%       |                   |               hd
%       _____________________               -


% Not shown in diagram:
% t_ft - axial thickness of the float top plate
% t_fb - axial thickness of the float bottom plate
% t_fc - circumferential thickness of the float gussets
% t_fr - radial thickness of the float walls
% t_sr - radial thickness of the spar walls

%% Surface float
num_gussets = 24;

% float cross sectional and lateral area
A_f_c = pi*(D_f+D_s)*t_fr + num_gussets * t_fc * (D_f-D_s)/2; % checked
A_f_l = num_gussets * t_fr * t_f;

% Volume of the surface float: (for material purposes)
V_top_plate = pi * (D_f/2)^2 * t_ft;
V_bot_plate = pi * (D_f/2)^2 * t_fb;
V_rims_gussets = A_f_c * h_f;
V_sf_m = V_top_plate + V_bot_plate + V_rims_gussets;

m_f_m = V_sf_m * rho_m(M); % mass of float material without ballast

% Volume of the surface float: (for displacement purposes)
A_f = pi/4 * (D_f^2 - D_s^2);
V_f_d = A_f * T_f;
m_f_tot = V_f_d * rho_w;
m_f_b = m_f_tot - m_f_m; % mass of ballast on float
V_f_b = m_f_b / rho_b; % volume of ballast on float
V_f_tot = A_f * h_f; % total volume available on float
V_f_pct = V_f_b / V_f_tot; % percent of available volume used by ballast on float

%% Vertical column
h = D_s * 7; % scaling law to keep same proportions as RM3
draft_vc = h;

V_vc_d = pi * (D_s/2)^2 * h;            % volume for displacement purposes
D_ivc = D_s - 2*t_sr;                   % inner diameter
A_c_vc = pi * ((D_s/2)^2-(D_ivc/2)^2);  % cross sectional area
V_vc_m = A_c_vc * h;                    % volume for material purposes
I_vc = pi * (D_s^4 - D_ivc^4) / 64;     % area moment of inertia
A_l_vc = D_s * h;                       % lateral area

%% Reaction plate
draft_rp = t_r;
A_c_rp = pi * ((D_or/2)^2 - (D_s/2)^2); % cross sectional area
A_l_rp = pi * D_or * t_r;               % lateral area
V_rp = A_c_rp * t_r;                    % Volume (for displacement and material purposes)
I_rp = pi* (D_or^4)/64;

%% Totals
V_d = V_f_d + V_vc_d + V_rp; % Total Volume displaced
V_m = V_sf_m + V_vc_m + V_rp; % Total Volume of material
m_tot = V_m * rho_m(M);       % Total mass of material

A_c = [A_f_c, A_c_vc, A_c_rp];
A_lat_sub = [A_f_l A_l_vc A_l_rp];
r_over_t = [0,... % D_sft/(2*t_sf) 
            D_s/(2*t_sr),...
            0];%D_or/(2*t_r)];
I = [I_f, I_vc, I_rp];
draft = [draft_sf,draft_vc,draft_rp];

%% Buoyancy Calculations
Fb = rho_w * V_d * g;
Fg = m_tot * g;
B = Fb / Fg;

%% Metacentric Height Calculatons
KB = (t_f/2 + h + t_r)/2;	% center of buoyancy above the keel
KG =  t_f/2 + h + t_r;  	% center of gravity above the keel
I_f = pi/64 * D_f^4;
BM = I_f / V_d;            % V_d is the submerged/displaced volume
GM = KB + BM - KG;          % Metacentric Height


end

