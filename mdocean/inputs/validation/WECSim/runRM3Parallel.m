%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2014 National Renewable Energy Laboratory and National
% Technology & Engineering Solutions of Sandia, LLC (NTESS).
% Under the terms of Contract DE-NA0003525 with NTESS,
% the U.S. Government retains certain rights in this software.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a revised version of WEC-Sim/source/wecSimPCT.m that saves
% specific output variables with a given filename.
% fixme: this perhaps belongs better in userDefinedFunctions.m?

%% wecSimPCT
% WEC-Sim parallel computing toolbox executable
clear mcr imcr i j k l1 l2 m n name nseed kkk len numConditions
clear body waves simu output pto constraint ptoSim
global mcr

% open the local cluster profile
pool = parcluster('local');
totalNumOfWorkers=pool.NumWorkers;

% open the parallel pool, recording the time it takes
tic;
pool = gcp; % open the pool

fprintf('Opening the parallel pool took %g seconds.\n', toc)

if isempty(pool)
    error('Could not open parpool')
else
    fprintf("Number of workers: %i",pool.NumWorkers)
end

evalc('wecSimInputFile');

if isempty(simu.mcrMatFile) == 0
    load(simu.mcrMatFile);
else
    kkk=0;
    mcr.header = {'waves.height','waves.period'};
    if isempty(simu.mcrExcelFile) == 0
        mcr.waveSS = xlsread(simu.mcrExcelFile);
        mcr.waveSS(isnan(mcr.waveSS))=0;
        for i=2:length(mcr.waveSS(:,1))
            for j=2:length(mcr.waveSS(1,:))
                if (mcr.waveSS(i,j)>0)
                    kkk = kkk+1;
                    mcr.cases(kkk,1) = mcr.waveSS(i,1);
                    mcr.cases(kkk,2) = mcr.waveSS(1,j);
                end
            end
        end
    else
        for i=1:length(waves.height)
            for j=1:length(waves.period)
                kkk = kkk+1;
                mcr.cases(kkk,1) = waves.height(i);
                mcr.cases(kkk,2) = waves.period(j);
            end
        end
    end
    
    numConditions=2;
    if length(waves.phaseSeed)>1
        numConditions=numConditions+1;
        mcr.header{numConditions} = 'waves.phaseSeed';
        len = length(mcr.cases(:,1));
        for nseed=1:length(waves.phaseSeed)
            mcr.cases(len*(nseed-1)+1:len*(nseed-1)+len,1:numConditions-1) = mcr.cases(1:len,1:numConditions-1);
            mcr.cases(len*(nseed-1)+1:len*(nseed-1)+len,    numConditions) = waves.phaseSeed(nseed);
        end
    end
    
    if exist('pto','var')
        for n=1:size(pto,2)
            if (length(pto(n).damping)>1 || length(pto(n).stiffness)>1)
                numConditions=numConditions+2;
                name = sprintf('pto(%i).damping', n);
                mcr.header{numConditions-1} = name;
                name = sprintf('pto(%i).stiffness', n);
                mcr.header{numConditions  } = name;
                
                len = length(mcr.cases(:,1)); kkk = 0;
                for l2=1:length(pto(n).stiffness)
                    for l1=1:length(pto(n).damping)
                        kkk=kkk+1;
                        mcr.cases(len*(kkk-1)+1:len*(kkk-1)+len,1:numConditions-2) = mcr.cases(1:len,1:numConditions-2);
                        mcr.cases(len*(kkk-1)+1:len*(kkk-1)+len,  numConditions-1) = pto(n).damping(l1);
                        mcr.cases(len*(kkk-1)+1:len*(kkk-1)+len,    numConditions) = pto(n).stiffness(l2);
                    end
                end
            end
        end; clear i j k l1 l2 m n name nseed kkk len numConditions
    end
end

%% Execute wecSimPCT
pause(1)
delete savedLog*

% variables to save
STIFFNESS_COLUMN       = 4;
timesteps_per_period = mcr.cases(:,2) / simu.dt; 
P                      = zeros(length(mcr.cases(:,1)), 1);
force_pto              = zeros(length(mcr.cases(:,1)), 1);

float_amplitude        = zeros(length(mcr.cases(:,1)), 1);
spar_amplitude         = zeros(length(mcr.cases(:,1)), 1);

relative_amplitude     = zeros(length(mcr.cases(:,1)), 1);

float_amplitude_rms    = zeros(length(mcr.cases(:,1)), 1);
spar_amplitude_rms     = zeros(length(mcr.cases(:,1)), 1);
relative_amplitude_rms = zeros(length(mcr.cases(:,1)), 1);

float_amplitude_fund   = zeros(length(mcr.cases(:,1)), 1);
spar_amplitude_fund    = zeros(length(mcr.cases(:,1)), 1);
rel_amplitude_fund     = zeros(length(mcr.cases(:,1)), 1);

float_phase            = zeros(length(mcr.cases(:,1)), 1);
spar_phase             = zeros(length(mcr.cases(:,1)), 1);
rel_phase              = zeros(length(mcr.cases(:,1)), 1);

float_drag_force_fund  = zeros(length(mcr.cases(:,1)), 1);
spar_drag_force_fund   = zeros(length(mcr.cases(:,1)), 1);
float_drag_force_phase = zeros(length(mcr.cases(:,1)), 1);
spar_drag_force_phase  = zeros(length(mcr.cases(:,1)), 1);

float_pos_THD = NaN(length(mcr.cases(:,1)), 1);

% hardcoded corner sea states: [H, T] where H = Hs/sqrt(2)
H_all = mcr.cases(:,1);
T_all = mcr.cases(:,2);
corner_HT = [0.75/sqrt(2), 5.5;   % low Hs (0.75 m), low T
             0.75/sqrt(2), 11.5;  % low Hs (0.75 m), high T
             6.75/sqrt(2), 12.5;  % high Hs (6.75 m), low T
             6.75/sqrt(2), 15.5]; % high Hs (6.75 m), high T
% find matching indices in mcr.cases
corner_idx = zeros(size(corner_HT, 1), 1);
for ci = 1:size(corner_HT, 1)
    corner_idx(ci) = find(ismembertol(H_all, corner_HT(ci,1)) & ismembertol(T_all, corner_HT(ci,2)), 1);
end
corner_idx = corner_idx(:);
is_corner = false(length(mcr.cases(:,1)), 1);
is_corner(corner_idx) = true;

% preallocate cell arrays for timeseries at four corners
num_cases = length(mcr.cases(:,1));
float_accel_ts_cell = cell(num_cases, 1);
float_drag_ts_cell  = cell(num_cases, 1);
corner_N_per_T = round(timesteps_per_period(corner_idx));

parfor imcr=1:length(mcr.cases(:,1))
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    t = getCurrentTask();
    filename = sprintf('savedLog%03d.txt', t.ID);
    pctDir = sprintf('pctDir_%g', t.ID);
    getAttachedFilesFolder(pctDir);
    fileID = fopen(filename,'a');
    cleanupObj = onCleanup(@()cleanup_fcn(fileID,pctDir));
    fprintf(fileID,'wecSimPCT Case %g/%g on Worker Number %g/%g \n',imcr,length(mcr.cases(:,1)),t.ID,totalNumOfWorkers);
    % Run WEC-Sim
    try
        output = myWecSimFcn(imcr,mcr,pctDir,totalNumOfWorkers,X,p);   
    
        % extract signals over the last full period starting at a wave node,
        % so that the FFT phase is relative to the wave phase at t=0
        T_wave = mcr.cases(imcr,2);
        N_per_T = round(timesteps_per_period(imcr));
        wave_freq  = 2*pi/T_wave;
        k_last  = floor(simu.endTime / T_wave) - 1; % largest integer k such that k*T + T <= endTime
        i_start = round(k_last * T_wave / simu.dt) + 1; % 1-indexed; round() safe since T_wave/dt is always an integer
        i_end   = i_start + N_per_T - 1;
        power = output.ptos.powerInternalMechanics(i_start:i_end,3);
        F_PTO = output.ptos.forceInternalMechanics(i_start:i_end,3);
        float_pos = output.bodies(1).position(i_start:i_end,3);
        spar_pos  = output.bodies(2).position(i_start:i_end,3);
        rel_pos = float_pos - spar_pos;
        if isfield(output.bodies(1),'velocity') && isfield(output.bodies(2),'velocity')
            rel_vel = output.bodies(1).velocity(i_start:i_end,3) - ...
                      output.bodies(2).velocity(i_start:i_end,3);
        else
            % fallback when velocity is not logged; used only for sign inference
            rel_vel = gradient(rel_pos, simu.dt);
        end

        % include stiffness contribution from the separate spring in total PTO force
        F_spring = [];
        if isfield(output,'springs') && ~isempty(output.springs)
            if isfield(output.springs,'forceInternalMechanics')
                F_spring = output.springs.forceInternalMechanics(i_start:i_end,3);
            elseif isfield(output.springs,'forceTotal')
                F_spring = output.springs.forceTotal(i_start:i_end,3);
            end
        end
        if isempty(F_spring)
            K_spring = 0;
            if size(mcr.cases,2) >= STIFFNESS_COLUMN
                % mcr.cases(:,STIFFNESS_COLUMN) stores spring stiffness from mcr.header = ...,'stiffness'
                K_spring = mcr.cases(imcr,STIFFNESS_COLUMN);
            end
            force_sign = infer_force_sign(F_PTO, rel_vel);
            F_spring = force_sign * K_spring * rel_pos;
        end
        F_total = F_PTO + F_spring;
        F_drag_f = output.bodies(1).forceMorisonAndViscous(i_start:i_end,3);
        F_drag_s = output.bodies(2).forceMorisonAndViscous(i_start:i_end,3);

        % save specific output variables
        P(imcr) = mean(power);

        force_pto(imcr) = 1/2 * (max(F_total) - min(F_total));
        float_amplitude(imcr) = 1/2 * (max(float_pos) - min(float_pos));
        spar_amplitude(imcr)  = 1/2 * (max(spar_pos)  - min(spar_pos));
        relative_amplitude(imcr) = 1/2 * (max(rel_pos) - min(rel_pos));
        
        float_amplitude_rms(imcr) = rms( float_pos - mean(float_pos) );
        spar_amplitude_rms(imcr)  = rms( spar_pos  - mean(spar_pos) );
        relative_amplitude_rms(imcr) = rms( rel_pos - mean(rel_pos) );

        [float_amplitude_fund(imcr),...
         float_phase(imcr)] = get_fundamental(float_pos, wave_freq, simu.dt);
        [spar_amplitude_fund(imcr),...
         spar_phase(imcr)] = get_fundamental(spar_pos, wave_freq, simu.dt);
        [rel_amplitude_fund(imcr),...
         rel_phase(imcr)] = get_fundamental(rel_pos, wave_freq, simu.dt);

        [float_drag_force_fund(imcr), ...
         float_drag_force_phase(imcr)] = get_fundamental(F_drag_f, wave_freq, simu.dt);
        [spar_drag_force_fund(imcr), ...
         spar_drag_force_phase(imcr)]  = get_fundamental(F_drag_s, wave_freq, simu.dt);

        % compute float position THD for all sea states
        float_pos_THD(imcr) = compute_pos_thd(float_pos, wave_freq, simu.dt);

        % save timeseries for four corner sea states
        if is_corner(imcr)
            float_accel = output.bodies(1).acceleration(i_start:i_end,3);
            float_accel_ts_cell{imcr} = float_accel;
            float_drag_ts_cell{imcr} = F_drag_f;
        end
        
    catch ME
        warning(ME.identifier,'WecSim errored for sea state H=%.2f, T=%.1f: %s',...
            mcr.cases(imcr,1),mcr.cases(imcr,2),getReport(ME, 'extended', 'hyperlinks', 'off'));
        P(imcr) = NaN;
        force_pto(imcr) = NaN;

        float_amplitude(imcr) = NaN;
        spar_amplitude(imcr)  = NaN;
        relative_amplitude(imcr) = NaN;

        float_amplitude_rms(imcr) = NaN;
        spar_amplitude_rms(imcr)  = NaN;
        relative_amplitude_rms(imcr) = NaN;

        float_amplitude_fund(imcr) = NaN;
        spar_amplitude_fund(imcr) = NaN;
        rel_amplitude_fund(imcr) = NaN;

        float_phase(imcr) = NaN;
        spar_phase(imcr) = NaN;
        rel_phase(imcr) = NaN;

        float_drag_force_fund(imcr) = NaN;
        spar_drag_force_fund(imcr)  = NaN;
        float_drag_force_phase(imcr) = NaN;
        spar_drag_force_phase(imcr)  = NaN;
    end
    try; Simulink.sdi.clear; catch ME; warning(ME.identifier, 'Simulink.sdi.clear failed: %s', ME.message); end
end

B_p = mcr.cases(:,3);
K_p = mcr.cases(:,STIFFNESS_COLUMN);
dt_sim = simu.dt;
case_H = mcr.cases(:,1);
case_T = mcr.cases(:,2);

% assemble corner timeseries into matrices
max_N = max(corner_N_per_T);
num_corners = length(corner_idx);
float_accel_ts = zeros(max_N, num_corners);
float_drag_ts  = zeros(max_N, num_corners);
for ci = 1:num_corners
    ts_accel = float_accel_ts_cell{corner_idx(ci)};
    ts_drag  = float_drag_ts_cell{corner_idx(ci)};
    if ~isempty(ts_accel)
        float_accel_ts(1:length(ts_accel), ci) = ts_accel;
    end
    if ~isempty(ts_drag)
        float_drag_ts(1:length(ts_drag), ci) = ts_drag;
    end
end

var_names = wecsim_var_names();
save(output_filename, var_names{:})

clear imcr totalNumOfWorkers

%[DOCS AUTOGENERATED]
function cleanup_fcn(fileID,pctDir)
% Function cleanup_fcn
%
% :param fileID: fileID
% :param pctDir: pctDir
    fclose(fileID);
    try
        rmdir(pctDir, 's');
    end
end

function phase = get_phase(signal, N_per_T)
    [~,peak_idx] = max(signal);
    phase = (peak_idx - 1) / N_per_T * 2*pi;
end

function [fund,phase] = get_fundamental(signal,wave_freq,dt)
    Fs = 2*pi/dt;
    L = length(signal);
    n = L;
    Y = fft(signal,n);
    P2 = Y/L;
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    freqs = 0:(Fs/n):(Fs/2-Fs/n);
    idx_wave_freq = ismembertol(freqs, wave_freq);
    fund = abs(P1(idx_wave_freq));
    phase = angle(P1(idx_wave_freq));
end

function thd = compute_pos_thd(pos, wave_freq, dt)
    N = length(pos);
    Fs = 2*pi/dt;
    Y = fft(pos, N);
    P2 = Y/N;
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    freqs = (0:floor(N/2)) * (Fs/N);
    harmonic_numbers = round(freqs / wave_freq);
    amplitudes = abs(P1(1:length(freqs)));
    idx_fund = find(harmonic_numbers == 1, 1);
    if isempty(idx_fund) || amplitudes(idx_fund) == 0
        thd = NaN;
        return
    end
    fund = amplitudes(idx_fund);
    harmonics = amplitudes;
    harmonics(1) = 0;           % remove DC
    harmonics(idx_fund) = 0;    % remove fundamental
    thd = sqrt(sum(harmonics.^2)) / fund * 100;
end

function s = infer_force_sign(force_signal, relative_velocity)
    % Infer sign convention between PTO force and relative velocity.
    % Returns +1 for near-zero coupling so spring-force reconstruction falls
    % back to +K*relative_displacement when damping signal is negligible.
    rel_tol = 1e-12; % strict relative tolerance for sign-only inference
    denom = relative_velocity(:).' * relative_velocity(:);
    denom_tol = max(1, abs(denom)) * rel_tol;
    if denom <= denom_tol
        s = 1;
        return
    end
    proportionality = force_signal(:).' * relative_velocity(:) / denom;
    prop_tol = max(1, abs(proportionality)) * rel_tol;
    if abs(proportionality) <= prop_tol
        s = 1;
    else
        s = sign(proportionality);
    end
end
