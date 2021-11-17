
function LCOE = econ(p,m_tot, V_m, M, cost_m, N_WEC, P_elec)

structural_cost = [m_tot V_m m_tot] .* cost_m;
devicestructure =  N_WEC * structural_cost(M);

hr_per_yr = 8766;
aep = N_WEC * P_elec*hr_per_yr/1000; % annual energy production: W to kWh per year

capex = p.development + p.infrastructure + N_WEC*p.mooring + devicestructure + N_WEC*p.pto ...
        + p.profitmargin + N_WEC*p.installation + p.contingency; 
opex = N_WEC*p.operations + p.postinstall + p.shoreoperations + N_WEC*p.replacement ...
        + N_WEC*p.consumables + N_WEC*p.insurance;


LCOE = (p.FCR*capex + opex)/aep;

end