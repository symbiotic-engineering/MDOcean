
function [LCOE, Lt] = econ(m_tot, V_m, M, cost_m, N_WEC, i_PT, d_shore, FOS, P_elec)

structural_cost = [m_tot V_m m_tot] .* cost_m;
devicestructure =  N_WEC * structural_cost(M);

development     = 4553000;
infrastructure  = 990000;
mooring         = N_WEC * 525000;
pto             = N_WEC * 623000;
profitmargin    = 356000;
installation    = N_WEC * 5909000;
contingency     = 1590000;
operations      = N_WEC * 27000;
postinstall     = 710000;
shoreoperations = 142000;
replacement     = N_WEC * 54000;
consumables     = N_WEC * 8000;
insurance       = N_WEC * 227000;

hr_per_yr = 8766;
aep = N_WEC * P_elec*hr_per_yr/1000; % annual energy production: W to kWh per year

capex = development + infrastructure + mooring + devicestructure + pto ...
        + profitmargin + installation + contingency; 
opex = operations + postinstall + shoreoperations + replacement ...
        + consumables + insurance;
FCR=10.8;% 

LCOE = (capex*FCR + opex)/aep; % temporarily removed fcr

Lt = 20; % lifetime (years) RM3 ref

end