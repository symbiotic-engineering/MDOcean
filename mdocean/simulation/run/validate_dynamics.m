clear;close all;clc

%% Round 1: singlebody
b = var_bounds('wecsim');
X = [b.X_noms; 1];

% 1. drag off, wamit coeffs, wecsim geometry: should match wecsim very well <2%
p = parameters('wecsim');
p.C_d_float = 0;
p.C_d_spar = 0;
p.use_MEEM = false;
pct_error_baseline = run_dynamic_validation(X,p);

% 2. 1 but drag on: gives me % error that comes from drag
p = parameters('wecsim');
p.use_MEEM = false;
pct_error_drag = run_dynamic_validation(X,p);

% 3. drag back off but meem coffs: gives % error that comes from meem
p = parameters('wecsim');
p.C_d_float = 0;
p.C_d_spar = 0;
pct_error_meem_all = run_dynamic_validation(X,p);

%     3a. meem with N=50: gives meem truncation error
p = parameters('wecsim');
p.C_d_float = 0;
p.C_d_spar = 0;
p.harmonics = 50;
pct_error_meem_all = run_dynamic_validation(X,p);

%     3b. wamit with zero phase: gives meem phase error
% right now this requires manual fiddling with wamit coeffs

% 4. drag meem interaction
p = parameters('wecsim');
pct_error_total = run_dynamic_validation(X,p);

%% Round 2: multibody

% 1. drag off, wamit coeffs, wecsim geometry: should match wecsim very well <2%
% 2. 1 but drag on: gives me % error that comes from drag
% 3. 1 but meem coffs: gives % error that comes from meem
%     3a. meem truncation
%     3b. wamit with zero phase

%% Round 3: report geometry
p = parameters();
b = var_bounds();
X = [b.X_noms; 1];

% 1. drag on (perhaps sweep), meem coeffs, report geometry, try to match report 10%?

% bonus that I have't checked yet: irregular waves, force saturation


function pct_error = run_dynamic_validation(X,p)
    % override to have fewer sea states for the sake of fast debugging
%     p.Hs = [1 2];
%     p.T = [5 6];
%     p.JPD = [.25 .25; .25 25];

    wecsim_filename = run_wecsim_validation(p);
    pct_error = power_matrix_compare(X,p,wecsim_filename);
end

function output_filename = run_wecsim_validation(p)
    % p needs to be in the workspace for runRM3Parallel script to work right
    runRM3Parallel % this script uses p and modifies it, and saves output_filename to workspace

end