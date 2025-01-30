function [LCOE,capex,opex] = LCOE_from_capex_design_power(capex_design_dep, N_WEC, P_elec, FCR, efficiency)
    
    % costs taken from 'CBS (Total)' tab of the RM3 cost breakdown structure
    % https://catalog.data.gov/ne/dataset/reference-model-3-cost-breakdown-rm3-wave-point-absorber
    
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
    
    LCOE = (FCR*capex + opex)./AEP;

end