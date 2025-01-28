
function [LCOE, capex_design_dep, ...
            capex, opex, pto, devicestructure] = econ(mass_material, M, cost_perkg_mult, N_WEC, P_elec, FCR, ...
                                cost_perN_mult, cost_perW_mult, F_max, P_max, efficiency)

% costs taken from 'CBS (Total)' tab of the RM3 cost breakdown structure
% https://catalog.data.gov/ne/dataset/reference-model-3-cost-breakdown-rm3-wave-point-absorber
% with curve fits done in dev/design_cost_scaling.m

% structural cost per wec
alpha_struct = 0.481;
cost_per_kg = ( 1.64e6 + 1.31e6 * N_WEC^(-alpha_struct) ) / 687000 * cost_perkg_mult(M);
devicestructure = cost_per_kg * mass_material;

% PTO cost per wec
alpha_pto = 0.206;
pto_const =    92593 + 1051  *N_WEC^(-alpha_pto);
pto_power = ( 0.4454 + 0.9099*N_WEC^(-alpha_pto) ) * P_max * cost_perN_mult;
pto_force = ( 0.0648 + 0.0893*N_WEC^(-alpha_pto) ) * F_max * cost_perW_mult;
pto = pto_const + pto_power + pto_force;

% sum design-dependent cost per wec
capex_design_dep = devicestructure + pto; 

% design-independent cost per wec
% includes development, infrastructure, mooring, profitmargin,
% installation, contingency. (no decommissioning intentionally)
alpha_non_design = 0.741;
capex_non_design_dep = 12.68e6 * N_WEC^(-alpha_non_design) + 1.24e6;

% total capex
capex_per_wec = capex_design_dep + capex_non_design_dep;
capex = capex_per_wec * N_WEC;

% opex = operation, postinstall, replacement, consumables, and insurance
alpha_opex = 0.5567;
opex_per_wec = 1.193e6 * N_WEC^-alpha_opex;
opex = opex_per_wec * N_WEC;

% power and energy
hr_per_yr = 8766;
P_avg = N_WEC * P_elec * efficiency;
AEP = P_avg * hr_per_yr / 1000; % annual energy production: W to kWh per year

LCOE = (FCR*capex + opex)/AEP;

end