%% math to match the report while tweaking distance above the water
% above = 1.9534 : 0.00001 : 1.9535;
% T_f_2 = 5.2 - above;
% T_f_1 = 4 - above;
% D_f = 20;
% D_s = 6;
% A_f = pi/4 * (D_f^2 - D_s^2);
% V_f_cyl = A_f * T_f_1;                      % displaced volume of float: hollow cylinder portion
% V_f_fr = pi/12 * (T_f_2 - T_f_1) ...
%     * (D_f^2 + D_s^2 + D_f*D_s);            % displaced volume of float: non-hollow frustum portion
% V_f_fr_mid = pi/4 * D_s^2 * (T_f_2 - T_f_1);% displaced volume of float: center cylinder to subtract from frustum
% V_f_fru_hol = V_f_fr - V_f_fr_mid;          % displaced volume of float: hollow frustum portion
% V_f_d = V_f_cyl + V_f_fru_hol              % total displaced volume of float

%% math to match the WAMIT gdf
% above = 2;
% T_f_2 = 5 - above;
% T_f_1 = 4 - above;
% D_f = 20;
% D_s = 6;
% A_f = pi/4 * (D_f^2 - D_s^2);
% V_f_cyl = A_f * T_f_1;                      % displaced volume of float: hollow cylinder portion
% D_flat = 10;
% V_f_fr = pi/12 * (T_f_2 - T_f_1) ...
%     * (D_f^2 + D_flat^2 + D_f*D_flat);            % displaced volume of float: non-hollow frustum portion
% V_f_fr_mid = pi/4 * D_s^2 * (T_f_2 - T_f_1);% displaced volume of float: center cylinder to subtract from frustum
% V_f_fru_hol = V_f_fr - V_f_fr_mid;          % displaced volume of float: hollow frustum portion
% V_f_d = V_f_cyl + V_f_fru_hol              % total displaced volume of float

%% match to match the CAD
above = 2;
T_f_2 = 5.2 - above;
T_f_1 = 4 - above;
D_f = 20;
D_s = 6.5;
A_f = pi / 4 * (D_f^2 - D_s^2);
V_f_cyl = A_f * T_f_1;                      % displaced volume of float: hollow cylinder portion
D_flat = 6.5;
V_f_fr = pi / 12 * (T_f_2 - T_f_1) ...
    * (D_f^2 + D_flat^2 + D_f * D_flat);            % displaced volume of float: non-hollow frustum portion
V_f_fr_mid = pi / 4 * D_s^2 * (T_f_2 - T_f_1); % displaced volume of float: center cylinder to subtract from frustum
V_f_fru_hol = V_f_fr - V_f_fr_mid;          % displaced volume of float: hollow frustum portion
V_f_d = V_f_cyl + V_f_fru_hol;              % total displaced volume of float

%% decision: to match WAMIT, use the second chunk of code.
