function [p,T] = parameters(mode)

% mode = 'wecsim': use parameters corresponding to RM3.out in WEC-Sim
% mode = anything else or not provided: use parmaters corresponding to RM3 report (default)

if nargin<1
    mode = '';
end

ksi2pa = 6894757; % 1000 pounds/in2 to Pascal
in2m = 0.0254;  % inch to meter
yd2m = 0.9144;  % yard to meter
lb2kg = 1/2.2;  % pound to kilogram

if strcmpi(mode,'wecsim')
    T_s_over_D_s = 29/6;
    h_d_over_D_s = 0.1/6;
    T_f_2_over_h_f = 3/5;
    D_f_b_over_D_f = 10/20;
    T_f_1_over_T_f_2 = 2/3;
    D_f_in_over_D_s = 6/6;
else
    T_s_over_D_s = 35/6;
    h_d_over_D_s = 1*in2m/6;
    T_f_2_over_h_f = 3.2/5.2;
    D_f_b_over_D_f = 6.5/20;
    T_f_1_over_T_f_2 = 2/3.2;
    D_f_in_over_D_s = 6.5/6;
end

file = 'Humboldt_California_Wave Resource _SAM CSV.csv';
jpd = trim_jpd(readmatrix(file,'Range','A3'));
g = 9.8;
spar_exc = get_spar_exc(g);

cols  = {'name',  'name_pretty','value','subsystem','sweep',  'description','idx'};
types = {'string','string',     'cell', 'string',   'logical','string',     'cell'};
T = table('Size',[0 length(cols)],'VariableTypes',types);

% idx is for vector parameters only and determines which element of the
% vector should be used for normalization in the param sensitvities. It can
% be a cell containing 'min','max', or an integer index.

T = [T;
    ...% Environmental parameters
    table("rho_w","\rho_w",{1000},"site",false,"water density (kg/m3)",{''});
    ...%table("mu","\mu",1e-3,"site",false,"dynamic viscosity of water (Pa s)");
    table("g","g",{g},"site",false,"acceleration of gravity (m/s2)",{''});
    table("h","h",{100},"site",true,"water depth (m)",{''});
    table("JPD","JPD",{jpd(2:end,2:end)},"site",false,...
        "joint probability distribution of wave (%)",{''});
    table("Hs","H_s",{jpd(2:end,1)},"site",true,"significant wave height (m)",{'max'});
    table("T","T",{jpd(1,2:end)},"site",true,"wave energy period (s)",{'max'});
    table("Hs_struct","H_{s,struct}",{[5 7 9 11.22 9 7 5]*1.9*sqrt(2)},...
        "site",true,"100 year significant individual wave height (m)",{'max'});
    table("T_struct","T_{struct}",{[5.57 8.76 12.18 17.26 21.09 24.92 31.70]},...
        "site",true,"100 year wave peak period (s)",{'max'});
    ...% 100 year extreme heights and periods from Berg 2011 
    ...
    ...% Materials: [  Structural Steel ASTM-A36 (Ductile)
    ...%               Concrete (Brittle)
    ...%               Stainless Steel 316 (Ductile) ASTM-A240 ]
    table("sigma_y","\sigma_y",{[36,4.5,30]* ksi2pa},"structures",true,"yield strength (Pa)",{1});
    table("sigma_e","\sigma_e",{[58*.45, 0, 75*.45]*ksi2pa},"structures",true,"endurance limit (Pa)",{1});
    table("rho_m","\rho_m",{[7850 2400 7900]},"structures",true,"material density (kg/m3)",{1});
    table("E","E",{[200e9, 5000*sqrt(4.5*ksi2pa), 200e9]},"structures",true,...
            "young's modulus (Pa)",{1});
    table("cost_perkg_mult","cost_m",{[4.28, 125/yd2m^3/2400, 1.84/lb2kg]/4.28},"economics",true,"material cost ($/kg)",{1});
        ...% RM3 CBS sheet 1.4 average of cells F21, F34, F46, F56
        ...% https://www.concretenetwork.com/concrete-prices.html
        ...% https://agmetalminer.com/metal-prices/
    table("nu","\nu",{[0.36 0 0.29]},"structures",true,"Poisson's ratio (-)",{1});
    table("FOS_min","FOS_{min}",{1.5},"structures",true,"minimum FOS (-)",{''});
    ...
    ...% Thicknesses and structures: float
    table("t_f_t_over_t_f_b","t_{f,t}/t_{f,b}",{0.50/0.56},"structures",true,...
        "float top to bottom thickness ratio (-)",{''});
    table("t_f_r_over_t_f_b","t_{f,r}/t_{f,b}",{0.44/0.56},"structures",true,...
        "float radial to bottom thickness ratio (-)",{''});
    table("t_f_c_over_t_f_b","t_{f,c}/t_{f,b}",{0.44/0.56},"structures",true,...
        "float circumferential to bottom thickness ratio (-)",{''});
    table("D_f_tu","D_{f,tu}",{20 * in2m},"structures",true,...
        "float support tube diameter (m)",{''}); % 24 in p156 report, 20 in cad
    table("t_f_tu","t_{f,tu}",{.5 * in2m},"structures",true,...
        "float support tube thickness (m)",{''});
    table("w_over_h_stiff_f","w_{stiff,f}/h_{stiff,f}",{1/16},"structures",...
        true,"float stiffener width to height ratio (-)",{''});
    table("num_sections_f","N_{sect}",{12},"structures",false,"number of float sections (-)",{''});
    ...
    ... Thicknesses and structures: damping plate
    table("t_d_tu","t_{d,tu}",{1.00 * in2m},"structures",true,...
        "damping plate support tube radial wall thickness (m)",{''});
    table("D_d_tu","D_{d,tu}",{48.00 * in2m},"structures",true,...
        "damping plate support tube diameter (m)",{''});
    table("theta_d_tu","\theta_{d,tu}",{atan(17.5/15)},"structures",true,...
        "angle from horizontal of damping plate support tubes (rad)",{''});
    table("h_over_h1_stiff_d","h_{stiff,d}/h_{1,stiff,d}",{[12.5 0.5 22  1]/22},...
        "structures",true,"damping plate stiffener height (m)",{'max'});
    table("w_over_h1_stiff_d","w_{stiff,d}/h_{1,stiff,d}",{[.5 10 1 12]/22},...
        "structures",true,"damping plate stiffener width (m)",{'max'});
    table("FOS_mult_d","FOS_{mult,d}",{7.5},"structures",true,...
        "damping plate factor of safety multiplier (-)",{''});
    table("num_terms_plate","N_{plate}",{100},"structures",false,...
        "number of terms for damping plate concentrated load (-)",{''});
    table("radial_mesh_plate","N_{r,plate}",{20},"structures",false,...
        "number of radial mesh points for damping plate (-)",{''});
    table("num_stiff_d","N_{stiff,d}",{24},"structures",false,...
        "number of damping plate stiffeners (-)",{''});
    ...
    ...% Economics
    table("m_scale","m_{scale}",{1.1},"economics",false,... 
        "factor to account for mass of neglected stiffeners (-)",{''});
    table("power_scale","P_{scale}",{85.9/117.7},"economics",false,...
        "factor to scale power for validation tuning (-)",{''});
    table("FCR","FCR",{0.113},"economics",true,... 
        "fixed charge rate (-), see RM3 report p63",{''});
    table("N_WEC","N_{WEC}",{100},"economics",true,"number of WECs in array (-)",{''});
    table("LCOE_max","LCOE_{max}",{1},"economics",true,"maximum LCOE ($/kWh)",{''});
%     table("avg_power_min","P_{avg,elec,min}",{100},"economics",true,... 
%         "minimum average electrical power (W)",{''}); % set to a negative number (not zero) to disable constraint
    table("eff_array","\eta_{array}",{0.95*0.98},"economics",true,...
        "array availability and transmission efficiency (-)",{''});
    table("cost_perN_mult", "$/N",{1},"economics",true,...
        "cost per Newton multiplier (-), 0.1656 $/N from https://doi.org/10.1016/j.ifacol.2022.10.531",{''});
    table("cost_perW_mult", "$/W",{1},"economics",true,...
        "cost per Watt multiplier (-)",{''});
    ...
    ...% Geometric ratios of bulk dimensions
    table("D_d_min","D_{d,min}",{30},"geometry",true,... 
        "minimum damping plate diameter",{''});
    table("D_d_over_D_s","D_d/D_s",{30/6},"geometry",true,... 
        "normalized damping plate diameter (-)",{''});
    table("T_s_over_D_s","T_s/D_s",{T_s_over_D_s},"geometry",true,... 
        "normalized spar draft (-)",{''});
    table("h_d_over_D_s","h_d/D_s",{h_d_over_D_s},"geometry",true,... 
        "normalized damping plate thickness (-)",{''});
    table("T_f_2_over_h_f","T_{f,2}/h_f",{T_f_2_over_h_f},"geometry",true,... 
        "normalized float draft (-)",{''});
    table("T_f_1_over_T_f_2","T_{f,1}/T_{f,2}",{T_f_1_over_T_f_2},...
        "geometry",true,"normalized float draft slant (-)",{''});
    table("D_f_b_over_D_f","D_{f,b}/D_f",{D_f_b_over_D_f},"geometry",true,... 
        "normalized diameter of float bottom (-)",{''});
    table("D_f_in_over_D_s","D_{f,in}/D_s",{D_f_in_over_D_s},"geometry",...
        true,"ratio of float inner diameter to spar diameter (-)",{''});
    ...
    ...% Dynamics: device parameters
    table("C_d_float","C_{d,float}",{1},"dynamics",true,"coefficient of drag for float",{''});
    table("C_d_spar","C_{d,spar}",{1},"dynamics",true,"spar coefficient of drag",{''});
    table("eff_pto","\eta_{pto}",{0.80},"dynamics",true,"PTO efficiency (-)",{''});
    ...
    ...% Dynamics: simulation type
    table("control_type","control type",{'damping'},"dynamics",false,... 
        "reactive or damping",{''});
    table("use_MEEM","use_MEEM",{true},"dynamics",false,... 
        "whether to use MEEM for hydro coeffs (boolean)",{''});
    table("use_multibody","use_multibody",{false},"dynamics",false,... 
        "whether to use multibody dynamics (boolean)",{''});
    ...
    ...% Dynamics: numerics and convergence
    table("X_tol","X_{tol}",{1e-2},"dynamics",false,... 
        "max allowable iteration error on magnitude of amplitude (m)",{''});
    table("phase_X_tol","\angleX_{tol}",{deg2rad(3)},"dynamics",false,... 
        "max allowable iteration error on phase angle of amplitude (rad)",{''});
    table("max_drag_iters","max drag iters",{40},"dynamics",false,... 
        "max number of iterations for drag convergence (-)",{''});
    table("harmonics","harmonics",{10},"dynamics",false,... 
        "number of harmonics to use for MEEM (int)",{''});
    table("besseli_argmax","\mathrm{argmax(I}_\{nu}(z))",700.5,"dynamics",false,... 
        "max argument of besseli before Inf overflow occurs",{''});
    ...
    ...% Dynamics: hydro coefficient data for nominal design
    table("spar_excitation_coeffs","spar excitation coeffs",{spar_exc},...
        "dynamics",false,"spar excitation hydro coeffs from WAMIT for nominal RM3",{''});
    table("hydro","hydro",{readWAMIT(struct(),"rm3.out",[])},"dynamics",...
        false,"function from WECSim",{''});
    table("F_heave_mult","F_{heave,mult}",{0.98},"dynamics",true,... 
        "multiplier to make heave force match with validation (-)",{''})];
    % 1.925 is required to make WAMIT match tank test, and 0.857 is
    % required to make MEEM match WAMIT

T.Properties.VariableNames = cols;

p = cell2struct([T.value], T.name, 1);

if nargout > 1
    % if outputing the table (for sensitivities), compute numeric idx
    numeric_idx = ones(height(T),1);
    use_max = strcmp(T.idx,'max');
    use_min = strcmp(T.idx,'min');
    use_num = cellfun(@isnumeric,T.idx);
    [~,numeric_idx(use_max)] = cellfun(@max, T.value(use_max));
    [~,numeric_idx(use_min)] = cellfun(@min, T.value(use_min));
    numeric_idx(use_num) = [T.idx{use_num}];
    value_normalize = cellfun(@(x,i) x(numeric_idx(i)), T.value, num2cell(1:height(T)).','UniformOutput',false);

    T.index_normalize = numeric_idx;
    T.value_normalize = value_normalize;
end

end
