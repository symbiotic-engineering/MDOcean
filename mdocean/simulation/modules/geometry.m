function [V_d, m_m, m_f_tot, m_s_tot,...
         A_c, A_lat_sub, ...
         I, T, V_f_pct, V_s_pct, GM,...
         A_dt, L_dt, mass,...
         CB_f_from_waterline,CG_f_from_waterline] = geometry(D_s, D_f, D_f_in, D_f_b, ...
                                            T_f_1, T_f_2, h_f, h_s, h_fs_clear, D_f_tu, ...
                                            t_f_t, t_f_r, t_f_c, t_f_b, t_f_tu, t_s_r, t_d_tu, ...
                                            D_d, D_d_tu, theta_d_tu, T_s, h_d, t_d, ...
                                            h_stiff_f, w_stiff_f, num_sect_f, ...
                                            h_stiff_d, w_stiff_d, num_stiff_d, ...
                                            M, rho_m, rho_w, m_scale)

%% Variable Definitions
% D: diameter
% T: draft
% h: height
% t: material thickness
% _f: float
% _s: spar
% _d: damping plate
% _d_tu: damping plate tubular support
% _f_tu: float tubular support

%                /\
%               /  \ Lft
% h_fs_clear-  /    \                               -
%           | /  ___ \                              |
%           -/  |   | \                             |
%           /____Df____\                        -   |
%           |           |                       |   |
% ----------|           |---------  -   -       hf  |
%           |           |    Tf1    |   |       |   |
%            \         /      -    Tf2  |       |   |
%             \__Dfb__/             _   |       _   |
%               |   |                   |           |
%               |   |                   |           hs
%               |Ds |                   Ts          |
%               |   |                   |           |
%             / |   | \                 |           |
%        Ldt/   |   |   \               |           |
%         /     |   |     \             |           |
%       _________Dd__________           |   -       |
%       |                   |           |   hd      |
%       _____________________           -   -       -


% Not shown in diagram:
% D_f_in - inner diameter of float
% t_f_t - axial thickness of the float top plate
% t_f_b - axial thickness of the float bottom plate
% t_f_c - circumferential thickness of the float gussets
% t_f_r - radial thickness of the float walls
% t_s_r - radial thickness of the spar walls
% t_d_tu - radial thickness of damping plate support tube walls
% t_f_tu - radial thicknss of float support tube walls

%% Float
num_gussets = 2*num_sect_f;
num_gussets_loaded_lateral = 2;

% float cross sectional and lateral area for structural purposes
D_f_mean = (D_f + D_f_b)/2;
A_f_cross_top = pi * (D_f + D_f_in) * t_f_r + num_gussets * t_f_c * (D_f - D_f_in)/2; % a ring with diameter D_f, a ring with diameter D_f_in, and gussets
A_f_cross_bot = pi * (D_f_in)       * t_f_r + num_gussets * t_f_c * (D_f_mean - D_f_in)/2; % a ring with diameter D_s, and bottom part of gussets
A_f_l = num_gussets_loaded_lateral * t_f_c * T_f_2;

% float material volume and mass - see drawings on p168-169 of notebook 11/26/24
V_top_plate = pi * ( (D_f/2)^2 - (D_f_in/2)^2 ) * t_f_t;
V_bot_plate = pi * (D_f_b/2)^2 * t_f_b;
slant_height = sqrt((D_f/2 - D_f_b/2)^2 + (T_f_2 - T_f_1)^2); 
V_bot_slant = pi/2 * (D_f + D_f_b) * slant_height * t_f_b;       % lateral area of frustum
V_outer_rim = pi * D_f * t_f_r * (h_f - (T_f_2 - T_f_1));
V_inner_rim = pi * D_f_in * t_f_t * h_f;
A_gusset = (h_f - (T_f_2 - T_f_1)) * (D_f - D_f_in)/2 + (T_f_2 - T_f_1) * (D_f_mean - D_f_in)/2;
V_gussets = num_gussets * A_gusset * t_f_c;

% float support tubes for attaching PTO
num_float_tubes = 6;
horiz_leg_tube = D_f/2;
vert_leg_tube = (1 + D_f/(D_f-D_s)) * (h_fs_clear + h_s - h_f + T_f_2 - T_f_1 - T_s); % h_fs_clear is vertical clearance between float tubes and spar when at rest
L_ft = sqrt(horiz_leg_tube^2 + vert_leg_tube^2); % length of float tube
V_f_tubes = num_float_tubes * (D_f_tu^2 - (D_f_tu - 2*t_f_tu)^2) * pi/4 * L_ft;

% float stiffeners
num_stiff_per_sect_f = 2; % top and bottom - neglecting the side stiffeners
A_stiff_f = num_stiff_per_sect_f * num_sect_f * h_stiff_f * w_stiff_f;
len_stiff_f = (D_f - D_f_in)/2  * 0.75; % goes ~3/4 of the way along the float
V_stiff_f = len_stiff_f * A_stiff_f;

V_sf_m = V_top_plate + V_bot_plate + V_outer_rim + V_inner_rim ...
        + V_bot_slant  + V_gussets + V_f_tubes + V_stiff_f;

m_f_m = V_sf_m * rho_m(M) * m_scale;      % mass of float material without ballast

% float hydrostatic calculations
A_f = pi/4 * (D_f^2 - D_f_in^2);
V_f_cyl = A_f * T_f_1;                      % displaced volume of float: hollow cylinder portion
V_f_fr = pi/12 * (T_f_2 - T_f_1) ...
    * (D_f^2 + D_f_b^2 + D_f*D_f_b);        % displaced volume of float: non-hollow frustum portion
V_f_fr_mid = pi/4 * D_f_in^2 * (T_f_2 - T_f_1);% displaced volume of float: center cylinder to subtract from frustum
V_f_fru_hol = V_f_fr - V_f_fr_mid;          % displaced volume of float: hollow frustum portion
V_f_d = V_f_cyl + V_f_fru_hol;              % total displaced volume of float
m_f_tot = V_f_d * rho_w;

% ballast
m_f_b = m_f_tot - m_f_m;                % mass of ballast on float
V_f_b = m_f_b / rho_w;                  % volume of ballast on float
V_f_tot = V_f_d + A_f * (h_f - T_f_2);  % total volume available on float
V_f_pct = V_f_b / V_f_tot;              % percent of available volume used by ballast on float

I_f = pi/64 * D_f^4;                    % area moment of inertia of float

%% Spar (vertical column and damping plate)

V_vc_d = pi/4 * D_s^2 * (T_s-h_d);      % vertical column volume displaced
V_d_d = pi/4 * D_d^2 * h_d;             % damping plate volume displaced
V_s_d = V_vc_d + V_d_d;                 % spar volume displaced (both column and damping plate)
m_s_tot = rho_w * V_s_d;                % total spar mass

% vertical column material use
D_vc_i = D_s - 2 * t_s_r;                % spar column inner diameter
A_vc_c = pi/4 * (D_s^2 - D_vc_i^2);     % spar column cross sectional area
V_vc_m_body = A_vc_c * (h_s - h_d);
A_vc_caps = pi/4 * D_vc_i^2;
t_vc_caps = (1/2 + 2.5) * 0.0254; % middle is 2.5", top is 0.5"
V_vc_caps = A_vc_caps * t_vc_caps;
num_stiff_vc = 12;
A_vc_stiff = 0.658 + 2*0.652 + 0.350; % triangles
t_vc_stiff = 0.0254;
V_vc_m_stiff = num_stiff_vc * A_vc_stiff * t_vc_stiff;
V_vc_m = V_vc_m_body + V_vc_m_stiff + V_vc_caps;          % volume of column material

% damping plate material use
A_d = pi/4 * D_d^2;                     % damping plate itself
num_supports = 4;
L_dt = (D_d - D_s) / (2*cos(theta_d_tu));
D_dt_i = D_d_tu - 2 * t_d_tu;
A_dt = pi/4 * (D_d_tu^2 - D_dt_i^2);      % support tube area
num_unique_stiffeners = length(h_stiff_d)/2;
num_stiff_repeats = num_stiff_d / num_unique_stiffeners;
A_stiff_d = num_stiff_repeats * sum(h_stiff_d .* w_stiff_d);
V_d_m = A_d * t_d + num_supports * A_dt * L_dt + A_stiff_d * (D_d - D_s)/2;

% total spar material use and mass
m_vc_m = V_vc_m * rho_m(M) * m_scale;
m_d_m = V_d_m * rho_m(M) * m_scale;
m_s_m = m_vc_m + m_d_m;                 % mass of spar material

% spar ballast
m_s_b = m_s_tot - m_s_m;                % mass of spar ballast
V_s_b = m_s_b / rho_w;                  % volume of spar ballast
V_pto = .75 * pi * 3^2 * 12;            % volume of PTO (assume constant)
V_s_tot = pi/4 * D_s^2 * h_s - V_pto;   % total volume available on spar for ballast
V_s_pct = V_s_b / V_s_tot;              % percent of available volume used by ballast on spar            

I_vc = pi * (D_s^4 - D_vc_i^4) / 64;    % area moment of inertia
A_vc_l = 1/2 * pi * D_s * (T_s-h_d);    % lateral area

% Reaction plate
A_d_c = pi/4 * (D_d^2);         % cross sectional area
A_d_l = 1/2 * pi * D_d * h_d;           % lateral area
I_rp = pi * D_d^4 / 64;

%% Totals

A_c = [A_f_cross_top, A_vc_c, A_d_c];
A_lat_sub = [A_f_l A_vc_l A_d_l];
I = [I_f, I_vc, I_rp];
T = [T_f_2, T_s, h_d];          % drafts: used to calculated F_surge in dynamics.m
m_m = m_f_m + m_s_m;            % total mass of material

V_d = [V_f_d, V_vc_d, V_d_d];   % volume displaced
mass = [m_f_m, m_vc_m, m_d_m];  % material mass of each structure

%% Metacentric Height Calculation

% see dev/cob_com_frustum.mlx for derivation of float COB and COM
D_term_1 = -D_f_b^2 - 2*D_f_b*D_f + 3*D_f^2;
D_term_2 = 2*D_f^2 - 2*D_f_b^2;
D_term_3 = -6*D_f_in^2 + 3*D_f_b^2 + 2*D_f_b*D_f + D_f^2;
CB_f_integral = D_term_1*T_f_1^2 + D_term_2*T_f_1*T_f_2 + D_term_3*T_f_2^2;
CB_f_from_waterline = 1/V_f_d * pi/48 * CB_f_integral;

T_term_1 =    T_f_1^2 + 2*T_f_1*T_f_2 - 3*T_f_2^2;
T_term_2 =  2*T_f_1^2                 - 2 *T_f_2^2;
T_term_3 = -3*T_f_1^2 - 2*T_f_1*T_f_2 + 5 *T_f_2^2 - 12*T_f_2*h_f + 6*h_f^2;
T_term_4 =                                           12*T_f_2*h_f - 6*h_f^2;
T_term_5 = 4*(   T_f_1 -   T_f_2);
T_term_6 = 4*(-2*T_f_1 + 2*T_f_2 - 3*h_f);
T_term_7 = 12*h_f;
CG_f_num = T_term_1*D_f_b^2 + T_term_2*D_f_b*D_f + T_term_3*D_f^2 + T_term_4*D_f_in^2;
CG_f_den = T_term_5*D_f_b^2 + T_term_5*D_f_b*D_f + T_term_6*D_f^2 + T_term_7*D_f_in^2;
CG_f_from_waterline = CG_f_num / CG_f_den;

% centers of buoyancy, measured from keel (bottom of damping plate)
CB_f = T_s - CB_f_from_waterline;
CB_vc = h_d + (T_s-h_d)/2;
CB_d = h_d/2;
CBs = [CB_f, CB_vc, CB_d];

% centers of gravity, measured from keel - assume even mass distribution
CG_f = T_s - CG_f_from_waterline;
CG_vc = h_d + (h_s-h_d)/2;
CG_d = h_d/2;
CGs = [CG_f, CG_vc, CG_d];

KB = CBs * V_d' / sum(V_d);     % center of buoyancy above the keel
KG = CGs * mass' / sum(mass);   % center of gravity above the keel

BM = I_f / sum(V_d);            % moment due to buoyant rotational stiffness
GM = KB + BM - KG;              % metacentric height

end