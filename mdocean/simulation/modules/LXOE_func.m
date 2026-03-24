function [LXOE_matrix,LXOE] = LXOE_func(levelized_X_matrix,P_elec,efficiency)

levelized_X = mysum(levelized_X_matrix);
hr_per_yr = 8766;
P_weighted = P_elec .* p.JPD / 100 * in.eff_array;
AEP_matrix = P_weighted * in.N_WEC * hr_per_yr / 1000; % W to kWh per year, all wecs
P_avg = N_WEC * P_elec * efficiency;
AEP = P_avg * hr_per_yr / 1000; % annual energy production: W to kWh per year
LXOE = levelized_X ./ AEP;
LXOE_matrix = levelized_X_matrix ./ AEP_matrix; % from LCOE_from_capex_design_power.m: $/year / kWh/year = $/kWh