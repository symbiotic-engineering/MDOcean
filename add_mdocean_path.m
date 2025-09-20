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
cd(mdocean_folder)
clear path s MDOcean_folder mdocean_folder