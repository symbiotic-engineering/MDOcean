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

% External MATLAB dependencies (openflash, wec-sim, safe) are declared in
% mip.yaml and provided by mip: run commands with `mip project run '<target>'`,
% which loads them before running (see setup_mip). Because mdocean/ is added
% to the path here, after mip loads the packages, MDOcean's own code (e.g.
% run_MEEM, readWAMIT) takes precedence over same-named package functions.

parallel.defaultClusterProfile('local'); % use local instead of global to avoid interference between multiple CI runners

clear path s MDOcean_folder mdocean_folder
