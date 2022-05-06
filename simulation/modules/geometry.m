function [V_d, V_m, m_tot, m_float, h, A_c, ...
         A_lat_sub, r_over_t, I, draft, B_f, B_s, GM] = geometry(D_i, D_sft, ...
                                                  t_sft, t_sf, t_sfb, t_vc, ...
                                                  D_or, t_r, rho_m, M, rho_w)

% D_sft-diameter of the top surface float plate
% t_sft-thickness of the top surface float plate
% W_sf-width of the surface float column (rectangle)
% t_sf-thickness of the surface float column
% t_sfb-thickness of the bottom surface float plate

%% Surface float
t_f = D_sft / 4; % scaling law to keep same proportions as RM3
num_gussets = 24;

A_c_sf = pi*((D_sft/2)^2-(D_i/2)^2) + num_gussets * t_sf * (D_sft-D_i)/2;
A_l_sf = num_gussets * t_sf * t_f;

% Volume of the surface float: (for material purposes)
W_sf = pi * (D_sft + D_i); % circumference of inner + outer circle
V_top_plate = pi * (D_sft/2)^2 * t_sft;
V_bot_plate = pi * (D_sft/2)^2 * t_sfb;
V_gusset = num_gussets * t_sf * (D_sft-D_i)/2 * t_f;
V_rims = t_f * W_sf * t_sf;
V_sf_m = V_top_plate + V_rims + V_bot_plate + V_gusset;

m_float = V_sf_m * rho_m(M);
I_sf = pi/64 * D_sft^4;

% Volume of the surface float: (for displacement purposes)
V_sf_d = m_float / rho_w;
draft_sf = V_sf_d / (pi*((D_sft/2)^2 - (D_i/2)^2));
float_ratio_submerged = draft_sf / t_f;

%% Vertical column
h = D_i * 7; % scaling law to keep same proportions as RM3
draft_vc = h;

V_vc_d = pi * (D_i/2)^2 * h;            % volume for displacement purposes
D_ivc = D_i - 2*t_vc;                   % inner diameter
A_c_vc = pi * ((D_i/2)^2-(D_ivc/2)^2);  % cross sectional area
V_vc_m = A_c_vc * h;                    % volume for material purposes
I_vc = pi * (D_i^4 - D_ivc^4) / 64;     % area moment of inertia
A_l_vc = D_i * h;                       % lateral area

%% Reaction plate
draft_rp = t_r;
A_c_rp = pi * ((D_or/2)^2 - (D_i/2)^2); % cross sectional area
A_l_rp = pi * D_or * t_r;               % lateral area
V_rp = A_c_rp * t_r;                    % Volume (for displacement and material purposes)
I_rp = pi* (D_or^4)/64;

%% Totals
V_d = V_sf_d + V_vc_d + V_rp; % Total Volume displaced
V_m = V_sf_m + V_vc_m + V_rp; % Total Volume of material
m_tot = V_m * rho_m(M);       % Total mass of material

A_c = [A_c_sf, A_c_vc, A_c_rp];
A_lat_sub = [A_l_sf A_l_vc A_l_rp];
r_over_t = [0,... % D_sft/(2*t_sf) 
            D_i/(2*t_vc),...
            0];%D_or/(2*t_r)];
I = [I_sf, I_vc, I_rp];
draft = [draft_sf,draft_vc,draft_rp];

%% Buoyancy Calculations
Fb = rho_w * V_d * g;
Fg = m_tot * g;
B = Fb / Fg;

%% Metacentric Height Calculatons
KB = (t_f/2 + h + t_r)/2;	% center of buoyancy above the keel
KG =  t_f/2 + h + t_r;  	% center of gravity above the keel
BM = I_sf / V_d;            % V_d is the submerged/displaced volume
GM = KB + BM - KG;          % Metacentric Height


end

