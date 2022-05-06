
function [LCOE, scaledLCOE] = econ(m_tot, V_m, M, cost_m, N_WEC, P_elec, FCR)

structural_cost = [m_tot V_m m_tot] .* cost_m;
devicestructure =  N_WEC * structural_cost(M);

%Capital Cost: Device Costs
struc_assembly  = 4589696-devicestructure; 
pto             = 666180;
mooring         = 554872;

%Capital Cost: Balance of System Costs
development     = 1731579;
engimgmt        = 419469;
elecinfrastruct = 1830508;
plantcomiss     = 176324; 
siteaccess      = 121223;
asseminstall    = 1227951;

%Capital Cost: Financial Costs
contingency     = 551013;
insurconstruct  = 110203;
reserves        = 330608;

%O&M Costs
operations      = 888220;
maintenance     = 351125;

hr_per_yr = 8766;
aep = .269*(N_WEC * P_elec*hr_per_yr/1000); % annual energy production: W to kWh per year

capex = (struc_assembly + pto + mooring + development + engimgmt...
        + elecinfrastruct + plantcomiss + siteaccess + asseminstall...
        + contingency + insurconstruct + reserves);
    
scaledcapex = N_WEC*3951482.19 + 15935544.80;

opex = operations + maintenance;

scaledopex= N_WEC*42343.16 + 1197001; 

LCOE = (FCR*capex + opex)/aep;

scaledLCOE= (FCR*scaledcapex + scaledopex)/aep;
end