clear;close all;clc

filename = 'RM3-CBS.xlsx'; % spreadsheet containing RM3 "actual" power data

% inputs
p = parameters();
b = var_bounds();
X = [b.X_noms; 1];

% unsaturated power
actual_mech_unsat = readmatrix(filename,'Range','E73:S86','Sheet','Performance & Economics');
actual_elec_unsat = actual_mech_unsat * p.eff_pto;
[~, ~, P_matrix, ~, val] = simulation(X,p);
sim_elec_unsat = P_matrix/1000;
sim_mech_unsat = sim_elec_unsat / p.eff_pto;

% saturated power
v = validation_inputs();
p.power_max = v.power_max;
actual_elec_sat = readmatrix(filename,'Range','E97:S110','Sheet','Performance & Economics');
[~, ~, P_matrix] = simulation(X,p);
sim_elec_sat = P_matrix/1000;

% wecSim spar stationary
vars = {'P','float_amplitude','spar_amplitude'};
if p.C_d_float == 0
    spar_fixed = load('wecsim_sparfixed_floatcd0_ctrl-49aa381_4523abe',vars{:});
elseif p.C_d_float == 1
    spar_fixed = load('wecsim_sparfixed_floatcd1_92fac3d',vars{:});
else
    error('cant find wecsim data for this Cd')
end
spar_fixed_power = -reshape(spar_fixed.P, size(P_matrix)) / 1000;
spar_fixed_float_amplitude = reshape(spar_fixed.float_amplitude, size(P_matrix));

% wecSim spar moving
spar_moving = load('wecsim_power_sparfloatingcd5_floatcd0_multibody',vars{:});
spar_moving_power = -reshape(spar_moving.P, size(P_matrix)) / 1000;
spar_moving_power(spar_moving_power>1e6) = NaN;
spar_moving_float_amplitude = reshape(spar_moving.float_amplitude,size(P_matrix));

% wave resources
[T,H] = meshgrid(p.T, p.Hs);
JPD = p.JPD;
JPD_actual = readmatrix(filename,'Range','E24:S37','Sheet','Performance & Economics');
wave_resource_raw = 1030 * 9.8^2 / (64*pi) * T .* H.^2 / 1000;
wave_resource_sheet = readmatrix(filename,'Range','E49:S62','Sheet','Performance & Economics');
wave_resource_sheet(wave_resource_sheet == 0) = NaN;
wave_resource_sim = wave_resource_raw;
wave_resource_sim(JPD == 0) = NaN;

% variable manipulation before plotting
power_titles = {'RM3 Report','MDOcean',...
        'WecSim spar moving','WecSim spar fixed'};
    %'Actual: Saturated','Simulated: Saturated'};
pre_JPD_power        = actual_mech_unsat;
pre_JPD_power(:,:,2) = sim_mech_unsat;
%pre_JPD_power(:,:,3) = actual_elec_sat;
%pre_JPD_power(:,:,4) = sim_elec_sat;
pre_JPD_power(:,:,3) = spar_moving_power;
pre_JPD_power(:,:,4) = spar_fixed_power;

post_JPD_power = pre_JPD_power .* JPD / 100;

diameter = X(1);
CW = pre_JPD_power ./ wave_resource_raw;    % capture width
CW_max = 9.8 * T.^2 / (4*pi^2);             % max CW for linear hydrodynamics, axisymmetric body
CWR = CW / diameter;                        % capture width ratio
CW_to_CW_max = CW ./ CW_max;
CW_to_CW_max_zeroed = CW_to_CW_max;
JPD_rep = repmat(JPD,[1 1 4]);
CW_to_CW_max_zeroed(JPD_rep == 0) = NaN;

vars = pre_JPD_power;
%vars(:,:,:,2) = post_JPD_power;
vars(:,:,:,2) = CW;
% vars(:,:,:,4) = CWR;
vars(:,:,:,3) = CW_to_CW_max;
% vars(:,:,:,6) = CW_to_CW_max_zeroed;
vars(:,:,:,4) = pre_JPD_power ./ H.^2;

var_names = {'Unweighted Device Power Matrix (kW)',...
...%     'JPD-Weighted Power Matrix (kW)',...
     'Capture With (m)'...
...%     'Capture Width Ratio (-)',...
     'Capture Width / Max Capture Width (-)',...
...%     'Capture Width / Max Capture Width (-)',...
    'Unweighted Device Power Matrix per H^2 (kW/m^2)'};

% four figures showing weigted and unweighted power and CWR
for fig = 1:size(vars,4)
    var = vars(:,:,:,fig);
    figure
    for i = 1:4
        subplot(2,2,i)
        contourf(T,H,var(:,:,i))
        xlabel('T (s)')
        ylabel('Hs (m)')
        title(power_titles{i})
        colorbar
    end
    sgtitle(var_names{fig})

end

% resource figure
resource_titles = {'Raw Resource','Zeroed Resource (Actual)','Zeroed Resource (Sim)'};
wave_resources = wave_resource_raw;
wave_resources(:,:,2) = wave_resource_sheet;
wave_resources(:,:,3) = wave_resource_sim;
figure
for i = 1:3
    subplot(1,3,i)
    contourf(T,H,wave_resources(:,:,i))
    xlabel('Wave Period T (s)')
    ylabel('Wave Height Hs (m)')
    title(resource_titles{i})
    colorbar
end
sgtitle('Wave Resource (kW/m)')

% JPD figure
JPDs = JPD;
JPDs(:,:,2) = JPD_actual/100;
JPD_titles = {'Sim','Actual'};
figure
for i = 1:2
    subplot(1,2,i)
    contourf(T,H,JPDs(:,:,i))
    xlabel('Wave Period T (s)')
    ylabel('Wave Height Hs (m)')
    colorbar
    grid on
    title(JPD_titles{i})
end
sgtitle('JPD (-)')

% power error plot: vs RM3 report
figure                      % (sim - actual) / actual
pre_JPD_pct_error = 100 * (pre_JPD_power(:,:,2) - pre_JPD_power(:,:,1)) ./ pre_JPD_power(:,:,1);
%pre_JPD_pct_error(JPD==0) = NaN;
pre_JPD_pct_error(pre_JPD_pct_error==0) = NaN;
contourf(T,H,pre_JPD_pct_error,-10:10:100)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Unsaturated Power Percent Error to RM Report (%)')
colorbar
grid on

% power and amplitude error plot: vs WecSim
MDOcean_Wecsim_power_error = 100 * (pre_JPD_power(:,:,2) - pre_JPD_power(:,:,4)) ./ pre_JPD_power(:,:,4);
MDOcean_Wecsim_amplitude_error = 100 * (val.X - spar_fixed_float_amplitude) ./ spar_fixed_float_amplitude;
% set contour lines to make plot more readable
if p.C_d_float==0 && p.use_MEEM==false
    power_error_vals = 3;
    X_error_vals = 3;
elseif p.C_d_float==1 && p.use_MEEM==false
    power_error_vals = [-80 -50 -20 0 2:5 7 10 20];
    X_error_vals = [-50 -40 -20 -10 0:5 10];
elseif p.C_d_float==0 && p.use_MEEM==true
    power_error_vals = [-25:5:0 2 5];
    X_error_vals = -20:5:10;
elseif p.C_d_float==1 && p.use_MEEM==true
    power_error_vals = [-85 -50 -20 -14 -12 -10:5:20];
    X_error_vals = [-100 -10 -8 -7 -5 -2 0 2 5 10 20 30];
end
min_colorbar = min(min([MDOcean_Wecsim_power_error,MDOcean_Wecsim_amplitude_error],[],'all'),-5);
max_colorbar = max(max([MDOcean_Wecsim_power_error,MDOcean_Wecsim_amplitude_error],[],'all'),5);

figure
subplot 121
[c,h_fig] = contourf(T,H,MDOcean_Wecsim_power_error,power_error_vals);
clabel(c,h_fig);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Unsaturated Power')
colorbar
clim([min_colorbar max_colorbar])
colormap(bluewhitered)
grid on
subplot 122
[c,h_fig] = contourf(T,H,MDOcean_Wecsim_amplitude_error,X_error_vals);
clabel(c,h_fig);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Unsaturated Amplitude')
colorbar
clim([min_colorbar max_colorbar])
colormap(bluewhitered)
grid on
sgtitle('Percent Error to WecSim (%)')


% Bp and X (without power saturation)
if length(unique(val.B_p)) > 1
    figure
    contourf(T,H,val.B_p/1e6)
    xlabel('Wave Period T (s)')
    ylabel('Wave Height Hs (m)')
    title('B_p [MN/(m/s)]')
    colorbar
    grid on
end

figure
subplot 131
contourf(T,H,val.X)
title('MDOcean')
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
colorbar
grid on

subplot 132
contourf(T,H,spar_fixed_float_amplitude)
title('WEC-Sim spar fixed')
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
colorbar
grid on

subplot 133
contourf(T,H,spar_moving_float_amplitude)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('WEC-Sim spar moving')
colorbar
grid on

sgtitle('Float Amplitude (m)')

%% compare average power over all sea states in JPD
P_weighted_sim = sim_mech_unsat .* p.JPD / 100;
P_avg_sim = sum(P_weighted_sim(:))
P_weighted_actual = spar_fixed_power .* p.JPD / 100;
P_avg_actual = sum(P_weighted_actual(:))
pct_error_P_avg = 100 * (P_avg_sim - P_avg_actual) / P_avg_actual