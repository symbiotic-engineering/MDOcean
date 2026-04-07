function [LXOE_matrix,LXOE] = LXOE_func(levelized_X_matrix,P_elec,efficiency,JPD,N_WEC)

levelized_X = sum(levelized_X_matrix,'all','omitnan');
hr_per_yr = 8766;
P_weighted = P_elec .* JPD / 100 * efficiency;
AEP_matrix = P_weighted * N_WEC * hr_per_yr / 1000; % W to kWh per year, all wecs
%P_avg = N_WEC * P_elec * efficiency;
%AEP = P_avg * hr_per_yr / 1000; % annual energy production: W to kWh per year
AEP = sum(AEP_matrix,'all','omitnan');
LXOE = levelized_X ./ AEP;
LXOE_matrix = levelized_X_matrix ./ AEP_matrix; % from LCOE_from_capex_design_power.m: $/year / kWh/year = $/kWh