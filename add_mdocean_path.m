path = mfilename('fullpath'); % for regular scripts
if contains(path,'LiveEditorEvaluationHelper') % for matlab online
    path = matlab.desktop.editor.getActiveFilename;
    if isempty(path) % construct manually for jupyter notebooks (need to restart kernel if running this cell more than once)
        path = pwd; % assume running from tutorials folder
    end
end
s = split(which(path), filesep);
MDOcean_folder = strjoin(s(1:end-1), filesep);
mdocean_folder = [MDOcean_folder filesep 'mdocean'];
addpath(genpath(mdocean_folder))

rmpath(genpath([mdocean_folder '/simulation/modules/OpenFLASH'])) % prevent using OpenFLASH run_MEEM since it's not integrated yet

% allow WEC-Sim if it's installed in the parent directory of MDOcean or inside MDOcean
wecSim_folder_outside = [MDOcean_folder filesep '../WEC-Sim'];
wecSim_folder_inside = [MDOcean_folder filesep 'WEC-Sim'];
exist_outside = exist(wecSim_folder_outside,'dir');
exist_inside = exist(wecSim_folder_inside,'dir');
exist_vec = [exist_inside,exist_outside];
if any(exist_vec)
    folder_vec = {wecSim_folder_inside,wecSim_folder_outside};
    wecSim_folder = folder_vec{find(exist_vec,1)};
    wecSimSourceFolder = [wecSim_folder filesep 'source'];
    if isunix
        load_sl_glibc_patch % for linux, see https://www.mathworks.com/support/bugreports/2632298
    end
    set_param(0, 'ErrorIfLoadNewModel', 'off')
    addpath(genpath(wecSimSourceFolder))
    rmpath([wecSimSourceFolder '/functions/BEMIO/']) % prevent conflicts with MDOcean's modified readWAMIT
    clear wecSimSourceFolder wecSim_folder folder_vec
end

clear path s MDOcean_folder mdocean_folder wecSim_folder_outside wecSim_folder_inside exist_outside exist_inside exist_vec