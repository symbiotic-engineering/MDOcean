function [V_d, m_m, m_f_tot, ...
         A_c, A_lat_sub, r_over_t, ...
         I, T, V_f_pct, V_s_pct, GM, mass] = geometry(D_s, D_f, T_f, h_f, h_s, ...
                                            t_ft, t_fr, t_fc, t_fb, t_sr, t_dt, ...
                                            D_d, D_dt, theta_dt, T_s, h_d, ...
                                            M, rho_m, rho_w, m_scale)

%% Variable Definitions
% D: diameter
% T: draft
% h: height
% t: material thickness
% _f: float
% _s: spar
% _d: damping plate
% _dt: damping plate tubular support

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
%             / |   | \                 |           |
%        Ldt/   |   |   \               |           |
%         /     |   |     \             |           |
%       _________Dd__________           -   -       -
%       |                   |               hd
%       _____________________               -


% Not shown in diagram:
% t_ft - axial thickness of the float top plate
% t_fb - axial thickness of the float bottom plate
% t_fc - circumferential thickness of the float gussets
% t_fr - radial thickness of the float walls
% t_sr - radial thickness of the spar walls
% t_dt - radial thickness of damping plate support tube walls

%% Float
num_gussets = 24;
num_gussets_loaded_lateral = 2;

% float cross sectional and lateral area
A_f_c = pi * (D_f + D_s) * t_fr + num_gussets * t_fc * (D_f - D_s)/2;
A_f_l = num_gussets_loaded_lateral * t_fc * T_f;

% float material volume and mass
V_top_plate = pi * (D_f/2)^2 * t_ft;
V_bot_plate = pi * (D_f/2)^2 * t_fb;
V_rims_gussets = A_f_c * h_f;
V_sf_m = V_top_plate + V_bot_plate + V_rims_gussets;

m_f_m = V_sf_m * rho_m(M) * m_scale;      % mass of float material without ballast

% float hydrostatic calculations
A_f = pi/4 * (D_f^2 - D_s^2);
V_f_d = A_f * T_f;
m_f_tot = V_f_d * rho_w;
m_f_tot = 725833;
warning('Overriding mass properties')

% ballast
m_f_b = m_f_tot - m_f_m;        % mass of ballast on float
V_f_b = m_f_b / rho_w;          % volume of ballast on float
V_f_tot = A_f * h_f;            % total volume available on float
V_f_pct = V_f_b / V_f_tot;      % percent of available volume used by ballast on float

I_f = pi/64 * D_f^4;            % area moment of inertia of float

%% Spar (vertical column and damping plate)

V_vc_d = pi/4 * D_s^2 * T_s;    % vertical column volume displaced
V_d_d = pi/4 * D_d^2 * h_d;     % damping plate volume displaced
V_s_d = V_vc_d + V_d_d;         % spar volume displaced (both column and damping plate)
m_s_tot = rho_w * V_s_d;        % total spar mass
m_s_tot = 886687;

% vertical column material use
D_vc_i = D_s - 2 * t_sr;                 % spar column inner diameter
A_vc_c = pi/4 * (D_s^2 - D_vc_i^2);      % spar column cross sectional area
V_vc_m = A_vc_c * h_s;                   % volume of column material

% damping plate material use
A_d = pi/4 * D_d^2; % damping plate itself
num_supports = 4;
L_dt = D_d / (2*cos(theta_dt));
D_dt_i = D_dt - 2 * t_dt;
A_dt = pi/4 * (D_dt^2 - D_dt_i^2); % support tube area
V_d_m = A_d * h_d + num_supports * A_dt * L_dt;

% total spar material use and mass
m_vc_m = V_vc_m * rho_m(M) * m_scale;
m_d_m = V_d_m * rho_m(M) * m_scale;
m_s_m = m_vc_m + m_d_m;                 % mass of spar material

% spar ballast
m_s_b = m_s_tot - m_s_m;                % mass of spar ballast
V_s_b = m_s_b / rho_w;                  % volume of spar ballast
V_powertrain = .75 * pi * 3^2 * 12;
V_s_tot = pi/4 * D_s^2 * h_s - V_powertrain;      % total volume available on spar for ballast
V_s_pct = V_s_b / V_s_tot;              % percent of available volume used by ballast on spar            

I_vc = pi * (D_s^4 - D_vc_i^4) / 64;    % area moment of inertia
A_vc_l = 1/2 * pi * D_s * T_s;          % lateral area

% Reaction plate
A_d_c = pi/4 * (D_d^2 - D_s^2);         % cross sectional area
A_d_l = 1/2 * pi * D_d * h_d;           % lateral area
I_rp = pi * D_d^4 / 64;

%% Totals

A_c = [A_f_c, A_vc_c, A_d_c];
A_lat_sub = [A_f_l A_vc_l A_d_l];
r_over_t = [0,... % D_sft/(2*t_sf) 
            D_s/(2*t_sr),...
            0];%D_or/(2*t_r)];
I = [I_f, I_vc, I_rp];
T = [T_f, T_s, h_d];
m_m = m_f_m + m_s_m;                    % total mass of material

V_d = [V_f_d, V_vc_d, V_d_d]; % volume displaced
mass = [m_f_m, m_vc_m, m_d_m]; % material mass of each structure

%% Metacentric Height Calculation
% centers of buoyancy, measured from keel (bottom of damping plate)
CB_f = h_d + T_s - T_f/2;
CB_vc = h_d + T_s/2;
CB_d = h_d/2;
CBs = [CB_f, CB_vc, CB_d];

% centers of gravity, measured from keel (assume even mass distribution)
CG_f = h_d + T_s - T_f + h_f/2;
CG_vc = h_d + h_s/2;
CG_d = h_d/2;
CGs = [CG_f, CG_vc, CG_d];

KB = CBs * V_d' / sum(V_d);     % center of buoyancy above the keel
KG = CGs * mass' / sum(mass);   % center of gravity above the keel

BM = I_f / sum(V_d);            % moment due to buoyant rotational stiffness
GM = KB + BM - KG;              % metacentric height

end

