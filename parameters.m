function p = parameters()

ksi2pa = 6894757; % 1000 pounds/in2 to Pascal
in2m = 0.0254;  % inch to meter
yd2m = 0.9144;  % yard to meter
lb2kg = 1/2.2;  % pound to kilogram

% Materials: [  Structural Steel ASTM-A36 (Ductile)
%               Concrete (Brittle)
%               Stainless Steel 304 (Ductile) ]

file = 'C:\SAM\2020.11.29\wave_resource\Humboldt_California_Wave Resource _SAM CSV.csv';
jpd = readmatrix(file,'Range','A3');

p = struct( 'rho_w',    1000,...                % water density (kg/m3)
            'd_shore',  1000,...                % distance from shore (m)
            'g',        9.8,...                 % acceleration of gravity (m/s2)
            'JPD',      jpd(2:end,2:end),...    % joint probability distribution of wave (%)
            'Hs',       jpd(2:end,1),...        % wave height (m)
            'Hs_struct',10,...                  % 100 year wave height (m)
            'T',        jpd(1,2:end),...        % wave period (s)
            'T_struct', 20,...                  % 100 year wave period (s)
            'tfinal',   30,...                  % simulation duration (s)
            's0',       [0; 0],...              % initial state [m m/s]
            'dt',       0.01,...                % timestep (s)
            'sigma_y',  [36,4.5,30]* ksi2pa,... % yield strength (Pa)
            'rho_m',    [8000 2400 8000],...    % material density (kg/m3)
            'E',        [200e6, ...             % young's modulus (Pa)
                        5000*sqrt(4.5*ksi2pa),...
                        200e6],...              
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
            'i_PT',     1,...                   % powertrain index (-)
            'P_min',    70e3,...                % minimum power (W)
            'B_min',    1,...                   % minimum buoyancy ratio (-)
            'FOS_min',  3);                     % minimum FOS (-)
        
end
