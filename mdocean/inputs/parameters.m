function p = parameters()

ksi2pa = 6894757; % 1000 pounds/in2 to Pascal
in2m = 0.0254;  % inch to meter
yd2m = 0.9144;  % yard to meter
lb2kg = 1/2.2;  % pound to kilogram

% Materials: [  Structural Steel ASTM-A36 (Ductile)
%               Concrete (Brittle)
%               Stainless Steel 304 (Ductile) ]

file = 'Humboldt_California_Wave Resource _SAM CSV.csv';
jpd = trim_jpd(readmatrix(file,'Range','A3'));
g = 9.8;
spar_exc = get_spar_exc(g);

above = 2; % top of float above still water line

p = struct( 'rho_w',    1000,...                % water density (kg/m3)
            ...%'mu',       1e-3,...                % dynamic viscosity of water (Pa s)
            'g',        g,...                   % acceleration of gravity (m/s2)
            'h',        100,...                 % water depth (m)
            'JPD',      jpd(2:end,2:end),...    % joint probability distribution of wave (%)
            'Hs',       jpd(2:end,1),...        % wave height (m)
            'Hs_struct',11.9,...                % 100 year wave height (m)
            'T',        jpd(1,2:end),...        % wave period (s)
            'T_struct', 17.1,...                % 100 year wave period (s)
            'sigma_y',  [36,4.5,30]* ksi2pa,... % yield strength (Pa)
            'rho_m',    [8000 2400 8000],...    % material density (kg/m3)
            'E',        [200e9, ...             % young's modulus (Pa)
                        5000*sqrt(4.5*ksi2pa),...
                        200e9],...              
            'cost_m',   [4.28, ...        % material cost ($/kg)
                        125/yd2m^3/2400, ...
                        1.84/lb2kg],...
            ...% RM3 CBS sheet 1.4 average of cells F21, F34, F46, F56
            ...% https://www.concretenetwork.com/concrete-prices.html
            ...% https://agmetalminer.com/metal-prices/
            'm_scale',  1.25,...                % factor to account for mass of neglected stiffeners (-)
            't_ft',     0.50 * in2m,...         % float top thickness (m)
            't_fr',     0.44 * in2m,...         % float radial wall thickness (m)
            't_fc',     0.44 * in2m,...         % float circumferential gusset thickness (m)
            't_fb',     0.56 * in2m,...         % float bottom thickness (m)
            't_sr',     1.00 * in2m,...         % vertical column thickness (m)
            't_d_max',  1.00 * in2m,...         % max thickness of damping plate before making it hollow (m)
            't_dt',     1.00 * in2m,...         % damping plate support tube radial wall thickness (m)
            'D_dt',     48.00 * in2m,...        % damping plate support tube diameter (m)
            'theta_dt', atan(17.5/15),...       % angle from horizontal of damping plate support tubes (rad)
            'FOS_min',  1.5,...                 % minimum FOS (-)	
            'FCR',      0.113,...               % fixed charge rate (-), see RM3 report p63
            'N_WEC',    100,...                 % number of WECs in array (-)   
            'D_d_over_D_s', 30/6,...            % normalized damping plate diameter (-)
            'T_s_over_D_s', 29/6,... %35/6      % normalized spar draft (-)
            'h_d_over_D_s', 0.1/6,...%in2m/6    % normalized damping plate thickness (-)     
            'T_f_2_over_h_f', (5-above)/5,...   % normalized float draft (-)
            'T_f_1_over_T_f_2',(4-above)/(5-above),...% normalized float draft slant (-)
            'D_f_b_over_D_f', 10/20,...         % normalized diameter of float bottom (-)
            'C_d_float',0,...                   % coefficient of drag for float
            'C_d_spar', 5,...                   % spar coefficient of drag
            'LCOE_max', .5,...                  % maximum LCOE ($/kWh)
            'power_max', Inf,...                % maximum power (W)
            'eff_pto',   0.80,...               % PTO efficiency (-)
            'eff_array', 0.95*0.98,...          % array availability and transmission efficiency (-)
            'control_type', 'damping',...       % 'reactive' or 'constant impedance' or 'damping'
            'X_tol',     1e-2,...               % max allowable iteration error on magnitude of amplitude (m)
            'phase_X_tol',deg2rad(3),...        % max allowable iteration error on phase angle of amplitude (rad)
            'max_drag_iters',100,...            % max number of iterations for drag convergence (-)
            'use_MEEM',  false,...              % whether to use MEEM for hydro coeffs (boolean)
            'use_multibody', false,...          % whether to use multibody dynamics (boolean)
            'harmonics', 10,...                 % number of harmonics to use for MEEM (int)
            'spar_excitation_coeffs',spar_exc,...% spar excitation hydro coeffs from WAMIT for nominal RM3
            'hydro',readWAMIT(struct(),'rm3.out',[]) ); % function from WECSim
end
