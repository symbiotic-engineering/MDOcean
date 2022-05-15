function LCOE = econ(m_m, M, cost_m, N_WEC, P_elec, FCR)


structural_cost = m_m.* cost_m;
devicestructure = structural_cost(M);

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

capex = .32*(devicestructure + struc_assembly + pto + mooring + development + engimgmt...
        + elecinfrastruct + plantcomiss + siteaccess + asseminstall...
        + contingency + insurconstruct + reserves);
  %.32 is a scaling factor based on the cost of 1 WEC in SAM. We are testing different materials so we cannot simply use slope  
  
scaledcapex = N_WEC*capex+ 15935544.80;

%opex = operations + maintenance;

scaledopex= N_WEC*42343.16 + 1197001; 

%LCOE = (FCR*capex + opex)/aep;

LCOE= (FCR*scaledcapex + scaledopex)/aep;
end