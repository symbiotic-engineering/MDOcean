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
            'T',        6,...                   % wave period (s)
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

X = [ 20 10 30;     % outer diameter of float
      .3 .1 .5;     % inner diameter ratio of float
      30 15 45;     % outer diameter of reaction plate
      1 2 3;        % material
      10 1 20;      % Number of WECs in array
      1e7 5e6 5e7   % D_int
        ];
    
X_nom = X(:,1);
design_size = size(X);
var_num = design_size(1);
var_depth = design_size(2);
LCOE = X*inf;
opt_idx = zeros(var_num,1);
recommended = zeros(var_num,2);

number_runs = var_depth + (var_num-1)*(var_depth-1);
failed = cell(number_runs,1);
power = zeros(number_runs,1);
FOS = zeros(number_runs,1);
X_ins = zeros(number_runs, var_num);

design = 0;
for i = 1:var_num
    X_in = X_nom;
    for j = 1:var_depth
        if i == 1 || j~=1
            changed_entry = X(i,j);
            if ~isnan(changed_entry)
                design = design+1;
                X_in(i) = changed_entry;
                X_ins(design,:) = X_in;
                [LCOE_temp, ~, ~, B, FOS(design), power(design)] = simulation(X_in,p);
                [feasible, failed{design}] = is_feasible(power(design), B, FOS(design), p);
                if feasible
                    LCOE(i,j) = LCOE_temp;
                else
                    LCOE(i,j) = NaN;
                end
            end
        end
    end
    [~, opt_idx(i)] = min(LCOE(i,:));
    recommended(i,:) = [X(i,opt_idx(i)), opt_idx(i)];
end
[LCOE_op, ~, ~, B_op, FOS_op, power_op] = simulation(recommended(:,1),p);
[op_feasible, CC] = is_feasible(power_op, B_op, FOS_op, p);
% create table for display
var_names = {'D_sft',...    % outer diameter of float (m)
            'D_i/D_sft',... % inner diameter ratio of float (m)
            'D_or',...      % outer diameter of reaction plate (m)
            'M',...         % material (-)
            'N_WEC',...     % number of WECs in array (-)
            'D_int'};       % internal damping of controller (Ns/m)
results = array2table(X_ins, 'VariableNames', var_names);
LCOE = LCOE';
results = addvars(results, round(LCOE(LCOE~=Inf),1), round(power/1e3), FOS, failed, ...
    'NewVariableNames', {'LCOE ($/kWh)','Power (kW)','FOS (-)','ConstraintsFailed'});
disp(results)

function [LCOE, D_env, Lt, B, FOS, P_elec] = simulation(X, p)

% capital X is design variables in vector format (necessary for optimization)
% lowercase x is design variables in struct format (more readable)

x = struct( 'D_sft',X(1),...        % outer diameter of float (m)
            'D_i',  X(2)*X(1),...   % inner diameter of float (m)
            'L_sf', pi*X(1),...        % radial material thickness of float (m) 
            'D_or', X(3),...        % outer diameter of reaction plate (m)
            'M',    X(4),...        % material (-)
            'N_WEC',X(5),...        % number of WECs in array (-)
            'D_int',X(6));          % internal damping of controller (Ns/m)

[V_d, V_m, m_tot, m_float, h, t_f, A_c, A_lat_sub, r_over_t, I] = geometry(x.D_i, x.D_sft, p.t_sft, x.L_sf, ...
                         p.t_sf, p.t_sfb, p.t_vc, x.D_or, p.t_r, p.rho_m, x.M);
        
[F_hydro_heave, F_hydro_surge, F_ptrain, D_env, P_elec] = dynamicSimulation(x, p, m_float, t_f);

[B,FOS] = structures(V_d, m_tot, F_hydro_heave, F_hydro_surge, F_ptrain, x.M, h, p.rho_w, p.g, p.sigma_y, A_c, A_lat_sub, r_over_t, I, p.E);

[LCOE, Lt] = econ(m_tot, V_m, x.M, p.cost_m, x.N_WEC, p.i_PT, p.d_shore, FOS, P_elec);

end