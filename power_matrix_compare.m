clear;close all;clc

filename = 'RM3-CBS.xlsx'; % spreadsheet containing RM3 "actual" power data

% inputs
p = parameters();
b = var_bounds(p);
X = [b.X_noms; 1];

% unsaturated power
actual_mech_unsat = readmatrix(filename,'Range','E73:T87','Sheet','Performance & Economics');
actual_elec_unsat = actual_mech_unsat * p.pto_eff;
[~, P_var, ~, ~, ~, ~, ~, ~, P_elec, ~, P_matrix] = simulation(X,p);
sim_elec_unsat = P_matrix(:,5:end)/1000;

% saturated power
v = validation_inputs();
p.power_max = v.power_max;
actual_elec_sat = readmatrix(filename,'Range','E97:T111','Sheet','Performance & Economics');
[~, P_var, ~, ~, ~, ~, ~, ~, P_elec, ~, P_matrix] = simulation(X,p);
sim_elec_sat = P_matrix(:,5:end)/1000;

% wave resources
[T,H] = meshgrid(p.T(5:end), p.Hs);
JPD = p.JPD(:,5:end);
JPD_actual = readmatrix(filename,'Range','E24:T38','Sheet','Performance & Economics');
wave_resource_raw = 1030 * 9.8^2 / (64*pi) * T .* H.^2 / 1000;
wave_resource_sheet = readmatrix(filename,'Range','E49:T63','Sheet','Performance & Economics');
wave_resource_sim = wave_resource_raw;
wave_resource_sim(JPD == 0) = 0;

% variable manipulation before plotting
power_titles = {'Actual: Unsaturated','Simulated: Unsaturated','Actual: Saturated','Simulated: Saturated'};
pre_JPD_power        = actual_elec_unsat;
pre_JPD_power(:,:,2) = sim_elec_unsat;
pre_JPD_power(:,:,3) = actual_elec_sat;
pre_JPD_power(:,:,4) = sim_elec_sat;

post_JPD_power = pre_JPD_power .* JPD / 100;

CWR_unweighted = pre_JPD_power ./ wave_resource_raw;
CWR_weighted = post_JPD_power ./ wave_resource_raw;

vars = pre_JPD_power;
vars(:,:,:,2) = post_JPD_power;
vars(:,:,:,3) = CWR_unweighted;
vars(:,:,:,4) = CWR_weighted;

var_names = {'Unweighted Device Power Matrix (kW)','JPD-Weighted Power Matrix (kW)',...
    'Unweighted Capture Width Ratio (-)','JPD-Weighted Capture Width Ratio (-)'};

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
sgtitle('Wave Resource')

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
