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

load_sl_glibc_patch % for linux, see https://www.mathworks.com/support/bugreports/2632298

wecSimFolder = [MDOcean_folder filesep '../WEC-Sim'];
if exist(wecSimFolder,'dir')    
    wecSimSourceFolder = [wecSimFolder filesep 'source'];
    set_param(0, 'ErrorIfLoadNewModel', 'off')
    addpath(genpath(wecSimSourceFolder))
    rmpath([wecSimSourceFolder '/functions/BEMIO/readWAMIT.m'])
    clear wecSimSourceFolder
end

clear path s MDOcean_folder mdocean_folder wecSimFolder