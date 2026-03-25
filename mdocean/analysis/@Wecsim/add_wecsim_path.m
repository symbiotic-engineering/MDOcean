function add_wecsim_path()
% ADD_WECSIM_PATH Add WEC-Sim to the MATLAB path, including applying the Simulink glibc patch.
if isunix
    load_sl_glibc_patch % for linux, see https://www.mathworks.com/support/bugreports/2632298
    disp('Loaded glibc patch for Simulink')
end

s = split(mfilename('fullpath'), filesep);
MDOcean_folder = strjoin(s(1:end-4), filesep); % up from add_wecsim_path -> @Wecsim -> analysis -> mdocean -> MDOcean root
mdocean_folder = [MDOcean_folder filesep 'mdocean'];

% allow WEC-Sim if it's installed in the parent directory of MDOcean or inside MDOcean
wecSim_folder_outside = [MDOcean_folder filesep '../WEC-Sim'];
wecSim_folder_inside = [MDOcean_folder filesep 'WEC-Sim'];
exist_outside = exist(wecSim_folder_outside,'dir');
exist_inside = exist(wecSim_folder_inside,'dir');
exist_vec = [exist_outside, exist_inside];
if any(exist_vec)
    folder_vec = {wecSim_folder_outside, wecSim_folder_inside};
    wecSim_folder = folder_vec{find(exist_vec,1)};
    wecSimSourceFolder = [wecSim_folder filesep 'source'];
    set_param(0, 'ErrorIfLoadNewModel', 'off')
    addpath(genpath(wecSimSourceFolder))
    addpath([mdocean_folder filesep 'inputs' filesep 'validation' filesep 'WECSim'],'-begin') % make sure MDOcean's modified readWAMIT takes precedence over WecSim's
end
end
