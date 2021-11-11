
function LCOE = econ(p,m_tot, V_m, M, cost_m, N_WEC, P_elec)

structural_cost = [m_tot V_m m_tot] .* cost_m;
devicestructure =  N_WEC * structural_cost(M);

mooring         = N_WEC * 525000;
pto             = N_WEC * 623000;
installation    = N_WEC * 5909000;
operations      = N_WEC * 27000;
replacement     = N_WEC * 54000;
consumables     = N_WEC * 8000;
insurance       = N_WEC * 227000;

hr_per_yr = 8766;
aep = N_WEC * P_elec*hr_per_yr/1000; % annual energy production: W to kWh per year

capex = p.development + p.infrastructure + mooring + devicestructure + pto ...
        + p.profitmargin + installation + p.contingency; 
opex = operations + p.postinstall + p.shoreoperations + replacement ...
        + consumables + insurance;


LCOE = (p.FCR*capex + opex)/aep;

end