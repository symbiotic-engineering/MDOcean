Dd = 30;
g = 9.8;

% inputs: geometric ratios
D_over_Dd = 0.2;
Td_over_Dd = 1.2;
td_over_Dd = 0.00085;
rho_Dd3_over_m = 60;

% calculate Cm mass coefficient
r = D_over_Dd;
Ca = rho_Dd3_over_m * 1/3 * (.5*r^3 - r^2 + 1);
Cm = Ca + 1;

% calculate natural frequency wn
denom = Td_over_Dd + 1/r^2 * td_over_Dd;
wn2 = 1/denom * g / (Cm * Dd);
wn = sqrt(wn2);
T = 2*pi/wn
