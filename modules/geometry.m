function [V_d, V_m, m_tot, m_float, h, t_f, A_c, A_lat_sub, r_over_t, I, draft] = geometry(~, ~, ~, ...
                                ~, ~, ~, ~, ~, ~, ~)

%D_sft-diameter of the top surface float plate
%t_sft-thickness of the top surface float plate
%W_sf-width of the surface float column (rectangle)
%t_sf-thickness of the surface float column
%t_sfb-thickness of the bottom surface float plate
                  
%% Surface float
% D_sft=20;
% D_i=6;
% t_sf=0.011176;
% t_sft=0.0127;
% t_sfb=0.014224;
% t_vc=0.0254;
% D_or=30;
% t_r=0.0254;
% rho_m=7850;
% M=[1];

float_pct_submerged = 1/2; % guess for now
t_f = D_sft / 4; % scaling law to keep same proportions as RM3
draft_sf = float_pct_submerged * t_f;

A_c_sf = pi*((D_sft/2)^2-(D_i/2)^2) + 24 * t_sf * (D_sft-D_i)/2;
A_l_sf = 24 * t_sf * t_f;

%Volume of the surface float: (for displacement purposes)
V_sf_d = pi*((D_sft/2)^2 - (D_i/2)^2) * draft_sf;
%Volume of the surface float: (for material purposes)
W_sf = pi * (D_sft + D_i); % circumference of inner + outer circle
V_sf_m = (pi * (D_sft/2)^2 * t_sft) + t_f * W_sf * t_sf + (pi * (D_sft/2)^2 * t_sfb) + + 24 * t_sf * (D_sft-D_i)/2 * t_f;
V_sf_s=0.18*V_sf_m; %total volume of the supports in the float
V_sf_tot=V_sf_m+V_sf_s; %total volume of th float with supports accounted for
m_float = V_sf_tot * rho_m(M); %mass of the float, given material, M, density
I_sf = pi/64 * D_sft^4;

%% Vertical column
h = D_i * 7; % scaling law to keep same proportions as RM3
draft_vc = h;

V_vc_d = pi * (D_i/2)^2 * h;            % volume for displacement purposes
D_ivc = D_i - t_vc ;                  % inner diameter
A_c_vc = pi * ((D_i/2)^2-(D_ivc/2)^2);  % cross sectional area
V_vc_m = A_c_vc * (h+(0.045*h));         % volume for material purposes
I_vc = pi * (D_i^4 - D_ivc^4) / 64;     % area moment of inertia
A_l_vc = D_i * h;                       % lateral area
V_vc_stiff= 1.72*V_vc_m;                  %volume of vertical column stiffeners (ratio)
V_vc_tot=V_vc_m+V_vc_stiff;             %Total volume of vertical column for material purpose
%% Reaction plate
draft_rp = t_r;
A_c_rp = pi * ((D_or/2)^2 - (D_i/2)^2); % cross sectional area
A_l_rp = pi * D_or * t_r;               % lateral area
V_rp = A_c_rp * t_r;                    % Volume (for displacement and material purposes)
V_rp_st=0.81*V_rp;                      % Volume of the the reaction plate radial supports and tubular structure (ratio)
V_rp_tot=V_rp+V_rp_st;                  % total volume of the reaction plate 
I_rp=pi* ((D_or^4)/64);
%% Totals
V_d = V_sf_d + V_vc_d + V_rp; % Total Volume displaced
V_m = V_sf_tot + V_vc_tot + V_rp_tot; % Total Volume of material
m_tot = V_m * rho_m(M);       % Total mass of material

A_c = [A_c_sf, A_c_vc, A_c_rp];
A_lat_sub = [A_l_sf A_l_vc A_l_rp];
r_over_t = [0,... % D_sft/(2*t_sf) 
            D_i/(2*t_vc),...
            0];%D_or/(2*t_r)];
I = [I_sf, I_vc, I_rp];
draft = [draft_sf,draft_vc,draft_rp];

end

