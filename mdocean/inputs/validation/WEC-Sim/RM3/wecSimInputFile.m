% add the desired outputs to runRM3Parallel.m and then run runRM3Parallel while in
% the mdocean/inputs/validation/WEC-Sim/RM3 folder

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

[~, ~, P_matrix, ~, val] = simulation(X,p);
if regular
    [T,H] = meshgrid(p.T, p.Hs/sqrt(2)); % regular
else
    [T,H] = meshgrid(p.T/0.857, p.Hs); % irregular
end

mcr.header = {'waves.height','waves.period','pto(1).damping'};
idx = p.JPD ~= 0;
control = val.B_p(idx);
mcr.cases = [H(idx),T(idx),control];
save('mcrMDOcean.mat','mcr')

% filename to save
[~, status] = system('git status');
if ~contains(status,'working tree clean')
    error('you have  uncommitted changes, please commit so the wecsim settings can be referenced to the commit')
end
[~, git_output] = system('git rev-parse --short HEAD');
git_hash = git_output(1:end-1);
uuid = char(matlab.lang.internal.uuid());
C_d_s = num2str(p.C_d_spar);
C_d_f = num2str(p.C_d_float);
mb = num2str(p.use_multibody);

output_filename = ['wecsim_sparcd' C_d_s '_floatcd' C_d_f '_multibody_' mb ...
                    '_meem_' meem '_git_' git_hash '_uuid_' uuid];

%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
if p.use_multibody
    simu.simMechanicsFile = 'inputs/validation/WEC-Sim/RM3/RM3_translation.slx';      % Specify Simulink Model File
else
    simu.simMechanicsFile = 'inputs/validation/WEC-Sim/RM3/RM3_fixed.slx';
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
body(1) = bodyClass('inputs/validation/WEC-Sim/RM3/hydroData/rm3.h5');      
    % Create the body(1) Variable, Set Location of Hydrodynamic Data File 
    % and Body Number Within this File.   
body(1).geometryFile = 'inputs/validation/WEC-Sim/RM3/geometry/float.stl';    % Location of Geomtry File
body(1).mass = 'equilibrium';                   
    % Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    % Weight.
body(1).inertia = [20907301 21306090.66 37085481.11];  % Moment of Inertia [kg*m^2]     

% Spar/Plate
body(2) = bodyClass('inputs/validation/WEC-Sim/RM3/hydroData/rm3.h5'); 
body(2).geometryFile = 'inputs/validation/WEC-Sim/RM3/geometry/plate.stl'; 
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


body(1).quadDrag.cd = [0 0 p.C_d_float 0 0 0];
body(1).quadDrag.area = [0 0 pi*(10^2-3^2) 0 0 0];
body(2).quadDrag.cd = [0 0 p.C_d_spar 0 0 0];
body(2).quadDrag.area = [0 0 pi*15^2 0 0 0];


