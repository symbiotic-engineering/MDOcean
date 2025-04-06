%%%% PTO cost %%%%%%
% PTO cost scale with N_WEC (from curve fit)
% C_PTO_per_wec = C_PTO_const + C_PTO_coeff * N_WEC ^ -C_PTO_alpha
C_PTO_const = 0.280e6;
C_PTO_coeff = 0.344e6;
C_PTO_alpha = 0.206;

% cost that scales with force: generator + bearing + mounting
cost_force_N1  = 25740+88411+28600;
cost_force_N100 = 20037+66308+5720;

F_max_nom = 0.926e6;
C_force_coeff = (cost_force_N1 - cost_force_N100) / (F_max_nom * (1 - 100^(-C_PTO_alpha)))
C_force_const = (cost_force_N1 / F_max_nom) - C_force_coeff

% cost that scales with neither force nor power: riser cable and controller
cost_neither_N1  = 88000+5644;
cost_neither_N100 = 88000+5000;

C_neither_coeff = (cost_neither_N1 - cost_neither_N100) / ((1 - 100^(-C_PTO_alpha)))
C_neither_const = (cost_neither_N1) - C_neither_coeff

P_max_nom = 286000;
C_power_coeff = (C_PTO_coeff - C_neither_coeff - C_force_coeff*F_max_nom) / P_max_nom
C_power_const = (C_PTO_const - C_neither_const - C_force_const*F_max_nom) / P_max_nom


cost_force_pred = (C_force_const + C_force_coeff*[1 100].^-C_PTO_alpha) * F_max_nom;
cost_force_act = [cost_force_N1 cost_force_N100];

cost_neither_pred = C_neither_const + C_neither_coeff*[1 100].^-C_PTO_alpha;
cost_neither_act = [cost_neither_N1 cost_neither_N100];

cost_pto_pred = C_PTO_const + C_PTO_coeff*[1 100].^-C_PTO_alpha
cost_power_pred = cost_pto_pred - cost_neither_pred - cost_neither_pred

%%%%%% non design cost %%%%%%%%
clc
N_WEC = [1 10 50 100];
development     = [4553389 8773812 11003159 10820060] ./ N_WEC
infrastructure  = [990000 4860000 7566000 17310000] ./ N_WEC
mooring         = [524775 4722975 23614875 47229750] ./ N_WEC
profitmargin    = [356252 2561152 11323295 21921723] ./ N_WEC
installation    = [5908552 9081973 21531225 37859591] ./ N_WEC
decommissioning = 0;%installation; % decomissioning cost not used in LCOE
contingency     = [1589545 5561144 18827150 35435836] ./ N_WEC

non_design_capex = development + infrastructure + mooring + profitmargin + installation + decommissioning + contingency;


