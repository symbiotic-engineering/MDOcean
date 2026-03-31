% sweep_geoms  Run or reload the geometry sweep and generate all figures.
%
% This script is a thin wrapper around the SweepGeoms analysis class.
% To run the full simulation from scratch:
%   obj = SweepGeoms(); obj.run_all_from_analysis();
%
% To reload previously saved results and regenerate figures:
%   obj = SweepGeoms(); obj.run_all_from_load();

close all
% Navigate to the MDOcean repo root (required by GenericAnalysis)
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fullfile(script_dir, '..', '..');
cd(repo_root);
obj = SweepGeoms();
obj.run_all_from_analysis();