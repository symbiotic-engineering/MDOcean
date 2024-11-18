function [p,T] = parameters()

ksi2pa = 6894757; % 1000 pounds/in2 to Pascal
in2m = 0.0254;  % inch to meter
yd2m = 0.9144;  % yard to meter
lb2kg = 1/2.2;  % pound to kilogram

file = 'Humboldt_California_Wave Resource _SAM CSV.csv';
jpd = trim_jpd(readmatrix(file,'Range','A3'));
g = 9.8;
spar_exc = get_spar_exc(g);

above = 2; % top of float above still water line (m)

cols  = {'name',  'name_pretty','value','subsystem','sweep',  'description'};
types = {'string','string',     'cell', 'string',   'logical','string'};
T = table('Size',[0 length(cols)],'VariableTypes',types);

T = [T;
    ...% Environmental parameters
    table("rho_w","\rho_w",{1000},"site",false,"water density (kg/m3)");
    ...%table("mu","\mu",1e-3,"site",false,"dynamic viscosity of water (Pa s)");
    table("g","g",{g},"site",false,"acceleration of gravity (m/s2)");
    table("h","h",{100},"site",true,"water depth (m)");
    table("JPD","JPD",{jpd(2:end,2:end)},"site",false,...
        "joint probability distribution of wave (%)");
    table("Hs","H_s",{jpd(2:end,1)},"site",true,"wave height (m)");
    table("Hs_struct","H_{s,struct}",{11.9},"site",true,"100 year wave height (m)");
    table("T","T",{jpd(1,2:end)},"site",true,"wave period (s)");
    table("T_struct","T_{struct}",{17.1},"site",true,"100 year wave period (s)");
    ...
    ...% Materials: [  Structural Steel ASTM-A36 (Ductile)
    ...%               Concrete (Brittle)
    ...%               Stainless Steel 304 (Ductile) ]
    table("sigma_y","\sigma_y",{[36,4.5,30]* ksi2pa},"structures",true,"yield strength (Pa)");
    table("rho_m","\rho_m",{[8000 2400 8000]},"structures",true,"material density (kg/m3)");
    table("E","E",{[200e9, 5000*sqrt(4.5*ksi2pa), 200e9]},"structures",true,...
        "young's modulus (Pa)");
    table("cost_m","cost_m",{[4.28, 125/yd2m^3/2400, 1.84/lb2kg]},"economics",true,"material cost ($/kg)");
        ...% RM3 CBS sheet 1.4 average of cells F21, F34, F46, F56
        ...% https://www.concretenetwork.com/concrete-prices.html
        ...% https://agmetalminer.com/metal-prices/
    ...
    ...% Thicknesses and structures
    table("t_d_max","t_{d,max}",{1.00 * in2m},"structures",true,...
        "max thickness of damping plate before making it hollow (m)");
    table("D_dt","D_{dt}",{48.00 * in2m},"structures",true,...
        "damping plate support tube diameter (m)");
    table("theta_dt","\theta_{dt}",{atan(17.5/15)},"structures",true,...
        "angle from horizontal of damping plate support tubes (rad)");
    table("FOS_min","FOS_{min}",{1.5},"structures",true,"minimum FOS (-)");
    ...
    ...% Economics
    table("m_scale","m_{scale}",{1.25},"economics",false,...
        "factor to account for mass of neglected stiffeners (-)");
    table("FCR","FCR",{0.113},"economics",true,...
        "fixed charge rate (-), see RM3 report p63");
    table("N_WEC","N_{WEC}",{100},"economics",true,"number of WECs in array (-)");
    table("LCOE_max","LCOE_{max}",{.5},"economics",true,"maximum LCOE ($/kWh)");
    table("eff_array","\eta_{array}",{0.95*0.98},"economics",true,...
        "array availability and transmission efficiency (-)");
    ...
    ...% Geometric ratios
    table("D_d_min","D_{d,min}",{30},"geometry",true,...
        "minimum damping plate diameter");
    table("D_d_over_D_s","D_d/D_s",{30/6},"geometry",true,...
        "normalized damping plate diameter (-)");
    table("T_s_over_D_s","T_s/D_s",{29/6},"geometry",true,...
        "normalized spar draft (-)");
    table("h_d_over_D_s","h_d/D_s",{0.1/6},"geometry",true,...
        "normalized damping plate thickness (-)");
    table("T_f_2_over_h_f","T_{f,2}/h_f",{(5-above)/5},"geometry",true,...
        "normalized float draft (-)");
    table("T_f_1_over_T_f_2","T_{f,1}/T_{f,2}",{(4-above)/(5-above)},...
        "geometry",true,"normalized float draft slant (-)");
    table("D_f_b_over_D_f","D_{f,b}/D_f",{10/20},"geometry",true,...
        "normalized diameter of float bottom (-)");
    ...
    ...% Dynamics: device parameters
    table("C_d_float","C_{d,float}",{0},"dynamics",true,"coefficient of drag for float");
    table("C_d_spar","C_{d,spar}",{5},"dynamics",true,"spar coefficient of drag");
    table("power_max","power_{max}",{Inf},"dynamics",true,"maximum power (W)");
    table("eff_pto","\eta_{pto}",{0.80},"dynamics",true,"PTO efficiency (-)");
    ...
    ...% Dynamics: simulation type
    table("control_type","control type",{'damping'},"dynamics",false,...
        "reactive or constant impedance or damping");
    table("use_MEEM","use_MEEM",{true},"dynamics",false,...
        "whether to use MEEM for hydro coeffs (boolean)");
    table("use_multibody","use_multibody",{false},"dynamics",false,...
        "whether to use multibody dynamics (boolean)");
    ...
    ...% Dynamics: numerics and convergence
    table("X_tol","X_{tol}",{1e-2},"dynamics",false,...
        "max allowable iteration error on magnitude of amplitude (m)");
    table("phase_X_tol","\angleX_{tol}",{deg2rad(3)},"dynamics",false,...
        "max allowable iteration error on phase angle of amplitude (rad)");
    table("max_drag_iters","max drag iters",{100},"dynamics",false,...
        "max number of iterations for drag convergence (-)");
    table("harmonics","harmonics",{10},"dynamics",false,...
        "number of harmonics to use for MEEM (int)");
    table("besseli_argmax","\mathrm{argmax(I}_\{nu}(z))",700.5,"dynamics",false,...
        "max argument of besseli before Inf overflow occurs");
    ...
    ...% Dynamics: hydro coefficient data for nominal design
    table("spar_excitation_coeffs","spar excitation coeffs",{spar_exc},...
        "dynamics",false,"spar excitation hydro coeffs from WAMIT for nominal RM3");
    table("hydro","hydro",{readWAMIT(struct(),"rm3.out",[])},"dynamics",...
        false,"function from WECSim")];

T.Properties.VariableNames = cols;

p = cell2struct([T.value], T.name, 1);

end
