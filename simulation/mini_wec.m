%% Fill in these inputs
in2m = 0.0254;

x = struct( 'D_sft', 4 * in2m,...   % diameter of WEC (m)
            'D_int', 10,...         % damping of controller (Ns/m)
            'w_n',   6 );           % natural frequency of controller (rad/s)

Hs_min = 1 * in2m;                  % minimum wave height (m)
Hs_max = 6 * in2m;                  % maximum wave height (m)
T_min  = 1;                         % minimum wave period (s)
T_max  = 2;                         % maximum wave period (s)
F_max = 1e3;                        % maximum powertrain force (N)
m_float = .5;                       % mass of WEC (kg)   

%% No need to edit this code
n = 10;
Hs = linspace(Hs_min,Hs_max,n);
T  = linspace(T_min, T_max, n);
JPD = ones(n) / n^2;

p = struct( 'rho_w',    1000,...        % water density (kg/m3)
            'g',        9.8,...         % acceleration of gravity (m/s2)
            'JPD',      JPD,...         % joint probability distribution of wave (%)
            'Hs',       Hs,...          % wave height (m)
            'Hs_struct',Hs_max,...      % 100 year wave height (m)
            'T',        T,...           % wave period (s)
            'T_struct', T_max,...       % 100 year wave period (s)
            'F_max',    F_max);         % max powertrain force (N)

% draft and Vd are set to zero since we don't need to calculate structural forces
draft = 0;
V_d = 0;

[F_hydro_heave, ~, ...
	F_ptrain, ~, P_elec, P_matrix] = dynamics(x, p, m_float, V_d, draft)

[TT,HH] = meshgrid(T,Hs/in2m);
figure
contourf(TT,HH,P_matrix);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (in)')
title('Power (W)')
colorbar
