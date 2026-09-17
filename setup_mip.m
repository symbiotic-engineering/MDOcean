function setup_mip()
%SETUP_MIP  One-time setup of MDOcean's mip-managed MATLAB dependencies.
%
%   MDOcean declares its external MATLAB packages (openflash, wec-sim, safe)
%   in mip.yaml, with exact versions pinned in mip.lock (see
%   https://mip.sh/docs/projects). This script:
%
%     1. Installs the mip package manager (https://mip.sh) if missing.
%     2. Installs and loads the preview build of mip that supports projects
%        and lockfiles (mip-org/labs/mip@mep-0009). This step is temporary,
%        until MEPs 8 and 9 land in a mip release.
%     3. Builds the project environment ./.mip to match mip.lock exactly
%        ("mip project sync").
%
%   Run this once after cloning, and again whenever mip.yaml or mip.lock
%   change (rerunning is cheap and idempotent). After setup, run MATLAB
%   commands with the dependencies loaded via
%
%     mip project run '<script, function, or expression>'
%
%   which is how the calkit pipeline stages run. (Until MEPs 8/9 land in a
%   mip release, load the preview build first in each session:
%   `mip load mip-org/labs/mip`.)

if isempty(which('/mip')) % non-builtin function named mip
    install_mip();
end

% Always use the preview build that provides "mip project" (MEPs 8/9).
% TEMPORARY: once MEPs 8/9 land in a mip release, drop these two lines --
% the released mip will have "mip project" and this install/load goes away.
mip('install', 'mip-org/labs/mip@mep-0009');
mip('load', 'mip-org/labs/mip');
rehash; % make sure subsequent mip calls dispatch to the preview build

MDOcean_folder = fileparts(mfilename('fullpath'));
mip('project', 'sync', '--yes', '--directory', MDOcean_folder);

fprintf(['\nMDOcean dependencies are ready in %s\n' ...
         'Run commands with them loaded via: mip project run ''<target>''\n'], ...
        fullfile(MDOcean_folder, '.mip'));
end

function install_mip()
% Run the official installer, which is non-interactive under matlab -batch.
% Kept in its own function so that any early-exit `return` inside the
% installer script returns only from here, not from setup_mip.
eval(webread('https://mip.sh/install.txt'));
end
