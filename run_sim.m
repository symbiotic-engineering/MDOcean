clear;clc;close all

ksi2pa = 6894757; % 1000 pounds/in2 to Pascal
in2m = 0.0254;  % inch to meter
yd2m = 0.9144;  % yard to meter
lb2kg = 1/2.2;  % pound to kilogram

% Materials: [  Structural Steel ASTM-A36 (Ductile)
%               Concrete (Brittle)
%               Stainless Steel 304 (Ductile) ]

p = struct( 'rho_w',    1000,...                % water density (kg/m3)
            'd_shore',  1000,...                % distance from shore (m)
            'g',        9.8,...                 % acceleration of gravity (m/s2)
            'Hs',       1,...                   % wave height (m)
            'T',        5,...                   % wave period (s)
            'tfinal',   30,...                  % simulation duration (s)
            's0',       [0; 0],...              % initial state [m m/s]
            'dt',       0.01,...                % timestep (s)
            'sigma_y',  [36,4.5,30]* ksi2pa,... % yield strength (Pa)
            'rho_m',    [8000 2400 8000],...    % material density (kg/m3)
            'cost_m',   [0.86/lb2kg, ...        % material cost
                        125/yd2m^3, ...         % [$/kg $/m3 $/kg] 
                        1.84/lb2kg],...  
            ...% https://agmetalminer.com/metal-prices/
            ...% https://www.concretenetwork.com/concrete-prices.html
            't_sft',    0.50 * in2m,...         % float top thickness (m)
            't_sf',     0.44 * in2m,...         % float column thickness (m)
            't_sfb',    0.56 * in2m,...         % float bottom thickness (m)
            't_r',      1.00 * in2m,...         % reaction plate thickness (m)
            't_vc',     1.00 * in2m,...         % vertical column thickness (m)
            'd_WEC',    10,...                  % distance between WECs (m)
            'd_farm',   1000,...                % distance to farm (m)
            'i_PT',     1);                     % powertrain index (-)

X = [5 1 1*in2m 5 1 1 1000];

[LCOE, D_env, Lt, B, FOS] = simulation(X,p)

function [LCOE, D_env, Lt, B, FOS] = simulation(X, p)

% capital X is design variables in vector format (necessary for optimization)
% lowercase x is design variables in struct format (more readable)

x = struct( 'D_sft',X(1),...	% outer diameter of float (m)
            'D_i',  X(2),...    % inner diameter of float (m)
            'L_sf', X(3),...    % radial material thickness of float (m) 
            'D_or', X(4),...    % outer diameter of reaction plate (m)
            'M',    X(5),...    % material (-)
            'N_WEC',X(6),...    % number of WECs in array (-)
            'D_int',X(7));      % internal damping of controller (Ns/m)

[V_d, V_m, m_tot, m_float, h] = geometry(x.D_i, x.D_sft, p.t_sft, x.L_sf, ...
                         p.t_sf, p.t_sfb, p.t_vc, x.D_or, p.t_r, p.rho_m, x.M);
        
[F_max, D_env, P_elec] = dynamicSimulation(x, p, m_float);

[B,FOS] = structures(V_d, m_tot, F_max, x.M, x.D_i, x.D_sft, h, x.D_or, p.rho_w, p.g, p.sigma_y);

[LCOE, Lt] = econ(m_tot, V_m, x.M, p.cost_m, x.N_WEC, p.i_PT, p.d_shore, FOS, P_elec);

end
