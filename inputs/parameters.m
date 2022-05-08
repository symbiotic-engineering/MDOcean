function p = parameters()

ksi2pa = 6894757; % 1000 pounds/in2 to Pascal
in2m = 0.0254;  % inch to meter
yd2m = 0.9144;  % yard to meter
lb2kg = 1/2.2;  % pound to kilogram

% Materials: [  Structural Steel ASTM-A36 (Ductile)
%               Concrete (Brittle)
%               Stainless Steel 304 (Ductile) ]

file = 'Humboldt_California_Wave Resource _SAM CSV.csv';
jpd = readmatrix(file,'Range','A3');

p = struct( 'rho_w',    1000,...                % water density (kg/m3)
            'g',        9.8,...                 % acceleration of gravity (m/s2)
            'JPD',      jpd(2:end,2:end),...    % joint probability distribution of wave (%)
            'Hs',       jpd(2:end,1),...        % wave height (m)
            'Hs_struct',10,...                  % 100 year wave height (m)
            'T',        jpd(1,2:end),...        % wave period (s)
            'T_struct', 20,...                  % 100 year wave period (s)
            'sigma_y',  [36,4.5,30]* ksi2pa,... % yield strength (Pa)
            'rho_m',    [8000 2400 8000],...    % material density (kg/m3)
            'E',        [200e6, ...             % young's modulus (Pa)
                        5000*sqrt(4.5*ksi2pa),...
                        200e6],...              
            'cost_m',   [0.86/lb2kg, ...        % material cost
                        125/yd2m^3/2400, ...    % $/kg 
                        1.84/lb2kg],...  
            ...% https://agmetalminer.com/metal-prices/
            ...% https://www.concretenetwork.com/concrete-prices.html
            't_ft',     0.50 * in2m,...         % float top thickness (m)
            't_fr',     0.44 * in2m,...         % float radial wall thickness (m)
            't_fc',     0.44 * in2m,...         % float circumferential gusset thickness (m)
            't_fb',     0.56 * in2m,...         % float bottom thickness (m)
            't_sr',     1.00 * in2m,...         % vertical column thickness (m)
            'B_min',    1,...                   % minimum buoyancy ratio (-)
            'FOS_min',  1.5,...                 % minimum FOS (-)	
            'D_d_min',  30,...                  % minimum damping plate diameter
            'FCR',      0.108,...               % fixed charge rate (-)	
            'N_WEC',    100,...                 % number of WECs in array (-)   
            'D_d_over_D_s', 30/6,...            % normalized damping plate diameter (-)
            'T_s_over_D_s', 35/6,...            % normalized spar draft (-)
            'h_d_over_D_s', 1*in2m/6);          % normalized damping plate thickness (-)                           
        
end
