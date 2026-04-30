% add the desired outputs to runRM3Parallel.m and then run runRM3Parallel while in
% the mdocean/inputs/validation/WECSim/RM3 folder

regular = true;

if ~exist('p','var')
    warning('parameters were loaded from scratch')
    p = parameters();
end
if ~exist('X','var')
    warning('design variables were loaded from scratch')
    b = var_bounds();
    X = [b.X_noms; 1];
end

if ~p.use_MEEM
    meem = 'off';
else
    meem = num2str(p.harmonics);
end
p.use_MEEM = false; % use WAMIT coeffs to get control so it's truly optimal

[~, P_matrix, ~, val] = simulation(X,p);
if regular
    [T,H] = meshgrid(p.T, p.Hs/sqrt(2)); % regular
else
    [T,H] = meshgrid(p.T/0.857, p.Hs); % irregular
end

mcr.header = {'waves.height','waves.period','pto(1).damping','stiffness'};
idx = p.JPD > 0 & ~isnan(p.JPD);
control_B = val.B_p(idx);
control_K = -val.K_p(idx);
mcr.cases = [H(idx),T(idx),control_B,control_K];
mcr.jpd = reshape(p.JPD(idx), [], 1);
save('mcrMDOcean.mat','mcr')

% filename to save
[~, status] = system('git status');
if ~contains(status,'working tree clean')
    msg = ['you have  uncommitted changes, please commit so the wecsim ' ...
        'settings can be referenced to the commit. Status: ' status];
    err = MException('MDOcean:WecSim:uncommitted',msg);
    %throw(err)
end
[~, git_output] = system('git rev-parse --short HEAD');
git_hash = git_output(1:end-1);
C_d_s = num2str(p.C_d_spar);
C_d_f = num2str(p.C_d_float);
mb = num2str(p.use_multibody);
if p.T_s_over_D_s==29/6
    geom = 'wecsim';
elseif p.T_s_over_D_s==35/6
    geom = 'report';
else
    geom='custom';
end

output_filename = ['../results/Wecsim/wecsim_sparcd' C_d_s '_floatcd' C_d_f ...
                    '_multibody_' mb '_meem_' meem '_geom_' geom];

%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
if p.use_multibody
    simu.simMechanicsFile = 'inputs/validation/WECSim/RM3/RM3_translation.slx';      % Specify Simulink Model File
else
    simu.simMechanicsFile = 'inputs/validation/WECSim/RM3/RM3_fixed.slx';
end
simu.mode = 'accelerator';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer = 'off';                   % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 100;                    % Wave Ramp Time [s]
simu.endTime = 200;                     % Simulation End Time [s]
simu.solver = 'ode4';                   % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.01; 							% Simulation time-step [s]
simu.mcrMatFile = 'mcrMDOcean.mat';
simu.saveWorkspace = false;
simu.b2b = 1;

%% Wave Information 
% % noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       % Initialize Wave Class and Specify Type  

% % Regular Waves  
if regular
    waves = waveClass('regular');           % Initialize Wave Class and Specify Type  
else
    waves = waveClass('irregular');
    waves.spectrumType = 'PM';                % Specify Wave Spectrum Type
end
waves.waterDepth = p.h;
% waves.height = 2.5;                     % Wave Height [m]
% waves.period = 8;                       % Wave Period [s]

% % Regular Waves with CIC
% waves = waveClass('regularCIC');          % Initialize Wave Class and Specify Type                                 
% waves.height = 2.5;                       % Wave Height [m]
% waves.period = 8;                         % Wave Period [s]

% % Irregular Waves using PM Spectrum 
%  waves = waveClass('irregular');           % Initialize Wave Class and Specify Type
%  waves.height = 2.5;                       % Significant Wave Height [m]
%  waves.period = 8;                         % Peak Period [s]
 waves.direction=[0];

% % Irregular Waves using JS Spectrum with Equal Energy and Seeded Phase
% waves = waveClass('irregular');           % Initialize Wave Class and Specify Type
% waves.height = 2.5;                       % Significant Wave Height [m]
% waves.period = 8;                         % Peak Period [s]
% waves.spectrumType = 'JS';                % Specify Wave Spectrum Type
% waves.bem.option = 'EqualEnergy';         % Uses 'EqualEnergy' bins (default) 
% waves.phaseSeed = 1;                      % Phase is seeded so eta is the same

% % Irregular Waves using PM Spectrum with Traditional and State Space 
% waves = waveClass('irregular');           % Initialize Wave Class and Specify Type
% waves.height = 2.5;                       % Significant Wave Height [m]
% waves.period = 8;                         % Peak Period [s]
% waves.spectrumType = 'PM';                % Specify Wave Spectrum Type
simu.stateSpace = 1;                      % Turn on State Space
% waves.bem.option = 'Traditional';         % Uses 1000 frequnecies

% % Irregular Waves with imported spectrum
% waves = waveClass('spectrumImport');      % Create the Wave Variable and Specify Type
% waves.spectrumFile = 'spectrumData.mat';  % Name of User-Defined Spectrum File [:,2] = [f, Sf]

% % Waves with imported wave elevation time-history  
% waves = waveClass('elevationImport');          % Create the Wave Variable and Specify Type
% waves.elevationFile = 'elevationData.mat';     % Name of User-Defined Time-Series File [:,2] = [time, eta]

%% Body Data
% Float
body(1) = bodyClass('inputs/validation/WECSim/RM3/hydroData/rm3.h5');      
    % Create the body(1) Variable, Set Location of Hydrodynamic Data File 
    % and Body Number Within this File.   
body(1).geometryFile = 'inputs/validation/WECSim/RM3/geometry/float.stl';    % Location of Geomtry File
body(1).mass = 'equilibrium';                   
    % Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    % Weight.
body(1).inertia = [20907301 21306090.66 37085481.11];  % Moment of Inertia [kg*m^2]     

% Spar/Plate
body(2) = bodyClass('inputs/validation/WECSim/RM3/hydroData/rm3.h5'); 
body(2).geometryFile = 'inputs/validation/WECSim/RM3/geometry/plate.stl'; 
body(2).mass = 'equilibrium';                   
body(2).inertia = [94419614.57 94407091.24 28542224.82];

%% PTO and Constraint Parameters
% Floating (3DOF) Joint
constraint(1) = constraintClass('Constraint1'); % Initialize Constraint Class for Constraint1
constraint(1).location = [0 0 0];               % Constraint Location [m]

% Translational PTO
pto(1) = ptoClass('PTO1');                      % Initialize PTO Class for PTO1
pto(1).stiffness = 0;                           % PTO Stiffness [N/m]
% pto(1).damping = 1200000;                       % PTO Damping [N/(m/s)]
pto(1).location = [0 0 0];                      % PTO Location [m]

%% Drag (Morison)
D_s   = X(1);        % inner diameter of float (m)
D_f   = X(2);        % outer diameter of float (m)
T_f_2 = X(3);
D_f_in = p.D_f_in_over_D_s * D_s;
D_d = p.D_d_over_D_s * D_s;
T_s = p.T_s_over_D_s * D_s;
% body(1).quadDrag.cd = [0 0 p.C_d_float 0 0 0];
% body(1).quadDrag.area = [0 0 pi/4*(D_f^2-D_f_in^2) 0 0 0];
% body(2).quadDrag.cd = [0 0 p.C_d_spar 0 0 0];
% body(2).quadDrag.area = [0 0 pi/4*D_d^2 0 0 0];

max_w = 2*pi/min(p.T);
max_k = max_w^2/p.g;
min_wavelength = 2*pi/max_k;
max_element_length = min_wavelength/10;
numMorisonElementsFloat = ceil(D_f/max_element_length);
numMorisonElementsSpar  = ceil(D_d/max_element_length);

z_ME_float = -T_f_2 + val.CG_f;
z_ME_spar  = -T_s   + val.CG_s;

body(1).morisonElement = make_morison_struct(D_f, D_f_in, numMorisonElementsFloat, p.C_d_float, z_ME_float);
body(2).morisonElement = make_morison_struct(D_d, 0,      numMorisonElementsSpar,  p.C_d_spar,  z_ME_spar);

function s = make_morison_struct(D_out,D_in,n,C_d,z_pos)

    [x_pos,A_z]  = circle_strips(D_out,D_in,n);

    cd_mat          = repmat([0 0 C_d],[n,1]);
    unit_normal_vec = repmat([0 0 1],  [n,1]);
    y_z_pos         = repmat([0 z_pos],[n,1]);
    A_x_y           = repmat([0 0],    [n,1]);

    rgME = [x_pos.' y_z_pos];
    A    = [A_x_y   A_z.'];
    
    s = struct('option',1,'cd',cd_mat,'ca',zeros(n,3),...
               'area',A,'VME',zeros(n,1),'rgME',rgME,'z',unit_normal_vec);
end

function [x_vec,areas] = circle_strips(D_out,D_in,n)
    element_length = D_out / n;
    x_start = -D_out/2 + element_length/2;
    x_vec = x_start + element_length * (0:(n-1));

    x_starts = x_vec - element_length/2;
    x_ends   = x_vec + element_length/2;

    no_cutout    = x_starts >= D_in/2 | x_ends <= -D_in/2;
    all_cutout   = abs(x_starts) <  D_in/2 & abs(x_ends) <  D_in/2;
    left_cutout  = abs(x_starts) <  D_in/2 & abs(x_ends) >  D_in/2;
    right_cutout = abs(x_starts) >  D_in/2 & abs(x_ends) <  D_in/2;
    mid_cutout   = x_starts < -D_in/2 & x_ends > D_in/2;
    assert(all(no_cutout + all_cutout + left_cutout + right_cutout + mid_cutout == 1))

    areas = zeros(size(x_starts));
    areas(no_cutout)    = circle_strip_area_no_cutout( x_starts(no_cutout),  x_ends(no_cutout),  D_out/2);
    areas(all_cutout)   = circle_strip_area_all_cutout(x_starts(all_cutout), x_ends(all_cutout), D_out/2, D_in/2);
    areas(left_cutout)  = circle_strip_area_all_cutout(x_starts(left_cutout), D_in/2, D_out/2, D_in/2) + ...
                         circle_strip_area_no_cutout(D_in/2, x_ends(left_cutout), D_out/2);
    areas(right_cutout) = circle_strip_area_all_cutout(-D_in/2, x_ends(right_cutout), D_out/2, D_in/2) + ...
                          circle_strip_area_no_cutout(x_starts(right_cutout), -D_in/2, D_out/2);
    areas(mid_cutout)   = circle_strip_area_no_cutout(x_starts(mid_cutout), -D_in/2, D_out/2) + ...
                          circle_strip_area_all_cutout(-D_in/2, D_in/2, D_out/2, D_in/2) + ...
                          circle_strip_area_no_cutout(D_in/2, x_ends(mid_cutout), D_out/2);

    assert(ismembertol(sum(areas), pi/4 * (D_out^2-D_in^2)));
end

function areas = circle_strip_area_all_cutout(x_starts,x_ends,R_out,R_in)
    areas = circle_strip_area_no_cutout(x_starts,x_ends,R_out) - circle_strip_area_no_cutout(x_starts,x_ends,R_in);
end

function areas = circle_strip_area_no_cutout(x_starts,x_ends,R)

    same_sign = sign(x_starts) == sign(x_ends);

    x_larger = max(abs(x_starts),abs(x_ends));
    x_larger(x_larger>R) = R; % prevent tiny imaginary area that sometimes happens due to finite precision effects 
    x_smaller = min(abs(x_starts),abs(x_ends));

    % derivation at https://ocw.mit.edu/courses/18-01sc-single-variable-calculus-fall-2010/28c569c5d8e79b7e1a2be1755de42d25_MIT18_01SCF10_Ses70a.pdf
    area_x_to_center = @(x) R^2*asin(x/R) + x .* sqrt(R^2 - x.^2);

    areas_same_sign = area_x_to_center(x_larger) - area_x_to_center(x_smaller);
    areas_diff_sign = area_x_to_center(x_larger) + area_x_to_center(x_smaller);

    areas = zeros(size(x_starts));
    areas(same_sign)  = areas_same_sign(same_sign);
    areas(~same_sign) = areas_diff_sign(~same_sign);

end

