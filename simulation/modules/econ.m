
function [LCOE, capex, opex] = econ(m_m, M, cost_m, N_WEC, P_elec, FCR, cost_perN, F_max, efficiency)

structural_cost = m_m.* cost_m;
devicestructure =  N_WEC * structural_cost(M);

% costs taken from 'CBS (Total)' tab of the RM3 cost breakdown structure
% https://catalog.data.gov/ne/dataset/reference-model-3-cost-breakdown-rm3-wave-point-absorber
development     = 4553000;
infrastructure  = 990000;
mooring         = N_WEC * 525000;
pto             = N_WEC * (623000 + F_max * cost_perN);
profitmargin    = 356000;
installation    = 5909000;
contingency     = 1590000;

capex = development + infrastructure + mooring + devicestructure + pto ...
        + profitmargin + installation + contingency; 

operations      = N_WEC * 27000;
postinstall     = 710000;
shoreoperations = 142000;
replacement     = N_WEC * 54000;
consumables     = N_WEC * 8000;
insurance       = (.8 + .2*N_WEC) * 227000;

opex = operations + postinstall + shoreoperations + replacement ...
        + consumables + insurance;

hr_per_yr = 8766;
P_avg = N_WEC * P_elec * efficiency;
aep = P_avg * hr_per_yr / 1000; % annual energy production: W to kWh per year

LCOE = (FCR*capex + opex)/aep;

end