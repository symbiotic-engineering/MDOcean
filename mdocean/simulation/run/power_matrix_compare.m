function weighted_power_error = power_matrix_compare(X, p, wecsim_filename, report, override)

if nargin==0
    % inputs
    p = parameters();
    b = var_bounds();
    X = [b.X_noms; 1];

    if p.use_multibody
        wecsim_filename = 'wecsim_sparcd0_floatcd0_55e8584';
        %'wecsim_power_sparfloatingcd5_floatcd0_multibody'
    else % singlebody = spar fixed
        if p.C_d_float == 0
            wecsim_filename = 'wecsim_sparfixed_floatcd0_ctrl-49aa381_4523abe';
        elseif p.C_d_float == 1
            wecsim_filename = 'wecsim_sparfixed_floatcd1_92fac3d';
        else
            error('cant find wecsim data for this Cd')
        end
    end
end

if nargin<4
    report = false;
end

results_wecsim = load_wecsim_results(wecsim_filename, p);

results_mdocean = compute_mdocean_results(X,p);

if report
    results_RM3_report = load_RM3_report_results(p.eff_pto, override);
    results_actual = results_RM3_report;
    results_sim = results_wecsim;
    results_sim(2) = results_mdocean;
    actual_str = 'Actual: RM3 Report';
    sim_str = {'Sim: WecSim','Sim: MDOcean'};
else
    results_actual = results_wecsim;
    results_sim = results_mdocean;
    actual_str = 'Actual: WecSim';
    sim_str = {'Sim: MDOcean'};
end

%var_names = {'power_mech_unsat', 'power_elec_unsat', 'power_elec_sat', ...
%                 'T', 'H', 'JPD', 'float_amplitude', 'spar_amplitude', ...
%                 'relative_amplitude', 'PTO_damping','CW','CW_to_CW_max'};
vars_to_plot = {'power_mech_unsat','CW_to_CW_max','float_amplitude','relative_amplitude','spar_amplitude'};
comparison_plot(p.T, p.Hs, results_actual, results_sim, vars_to_plot, actual_str, sim_str)

% compare average power over all sea states in JPD
weighted_power_error = zeros([1,length(results_sim)]);
for i=1:length(results_sim)
weighted_power_error(i) = compute_weighted_percent_error(results_sim(i).power_mech_unsat, ...
                                                      results_actual.power_mech_unsat, p.JPD);
end

end

%%
function old_stuff()

var_names = {'Unweighted Device Power Matrix (kW)',...
...%     'JPD-Weighted Power Matrix (kW)',...
     'Capture With (m)'...
...%     'Capture Width Ratio (-)',...
     'Capture Width / Max Capture Width (-)',...
...%     'Capture Width / Max Capture Width (-)',...
    'Unweighted Device Power Matrix per H^2 (kW/m^2)'};

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


% % Bp (without power saturation)
% if length(unique(val.B_p)) > 1
%     figure
%     contourf(T,H,val.B_p/1e6)
%     xlabel('Wave Period T (s)')
%     ylabel('Wave Height Hs (m)')
%     title('B_p [MN/(m/s)]')
%     colorbar
%     grid on
% end

end

function weighted_error = compute_weighted_percent_error(sim, actual, weights)
    sim_weighted = sim .* weights / 100;
    sim_avg = sum(sim_weighted(:));

    actual_weighted = actual .* weights / 100;
    actual_avg = sum(actual_weighted(:));
    weighted_error = compute_percent_error_matrix(actual_avg, sim_avg);
end

function comparison_plot(T, H, actual, sim, vars_to_plot, actual_str, sim_str)
    num_subplots = 2*length(sim)+1;
    for var_idx = 1:length(vars_to_plot)
        var_name = vars_to_plot{var_idx};
        figure
        % actual
        subplot(1,num_subplots,1)
        contour_plot(T,H, actual.(var_name), actual_str);

        for i=1:length(sim)
            % sim
            subplot(1,num_subplots,1+i)
            contour_plot(T,H,sim(i).(var_name), sim_str(i));
            % error
            error = compute_percent_error_matrix(actual.(var_name), sim(i).(var_name));
            subplot(1,num_subplots,1+length(sim)+i)
            error_levels = 10;
            error_plot(T,H,error,'Percent Error',error_levels);
        end
        
        title_split = strsplit(var_name,'_');
        title_str = cellfun(@mlreportgen.utils.capitalizeFirstChar,title_split,'UniformOutput',false);
        title_str_spaces = strcat(title_str," ");
        sgtitle(horzcat(title_str_spaces{:}))
    end
end

function error_plot(T, H, error, error_title, error_values)
    [c,h_fig] = contour_plot(T, H, error, error_title, error_values);
    clabel(c,h_fig);
    %clim([min_colorbar max_colorbar])
    if ~isempty(h_fig)
        colormap(h_fig.Parent,bluewhitered)
    end
end

function [c,h_fig] = contour_plot(T, H, Z, Z_title, Z_levels)
    if nargin<5
        Z_levels = 10;
    end

    if any(isfinite(Z),'all')
        [c,h_fig] = contourf(T,H,Z,Z_levels);
    else
        c = []; h_fig = [];
    end
    title(Z_title)
    xlabel('Wave Period T (s)')
    ylabel('Wave Height Hs (m)')
    colorbar
    grid on
end

function results = load_RM3_report_results(eff_pto, override)

    report_filename = 'RM3-CBS.xlsx'; % spreadsheet containing RM3 "actual" power data
    sheet = 'Performance & Economics';

    power_mech_unsat = readmatrix(report_filename,'Range','E73:S86',...
                                    'Sheet',sheet);
    
    
    Hs = readmatrix(report_filename,'Range','D73:D86','Sheet',sheet);
    Te = readmatrix(report_filename,'Range','E72:S72','Sheet',sheet);

    power_elec_sat = readmatrix(report_filename,'Range','E97:S110',...
                                    'Sheet',sheet);

    JPD = readmatrix(report_filename,'Range','E24:S37','Sheet',sheet)/100;

    wave_resource_sheet = readmatrix(report_filename,'Range','E49:S62','Sheet',sheet);
    wave_resource_sheet(wave_resource_sheet == 0) = NaN;

    if override
        power_mech_unsat = power_mech_unsat(1:2,1:2);
        Hs = Hs(1:2); Te = Te(1:2);
        JPD = JPD(1:2,1:2);
        wave_resource_sheet = wave_resource_sheet(1:2,1:2);
    end

    power_elec_unsat = power_mech_unsat * eff_pto;
    [T,H] = meshgrid(Te,Hs);

    results = assemble_results_struct(size(power_mech_unsat),...
                                        'power_mech_unsat',power_mech_unsat, ...
                                        'power_elec_unsat',power_elec_unsat, ...
                                        'power_elec_sat',power_elec_sat, ...
                                        'T',T, 'H',H, 'JPD',JPD);
end

function results = compute_mdocean_results(X,p)

    % unsaturated power
    [~, ~, P_matrix, ~, val] = simulation(X,p);
    power_elec_unsat = P_matrix/1000;
    power_mech_unsat = power_elec_unsat / p.eff_pto;
    
    % saturated power
    [~, ~, P_matrix] = simulation(X,p);
    power_elec_sat = P_matrix/1000;

    results = assemble_results_struct(size(power_mech_unsat),...
                                      'power_mech_unsat', power_mech_unsat, ...
                                      'power_elec_unsat',power_elec_unsat, ...
                                      'power_elec_sat', power_elec_sat, ...
                                      'T',p.T, 'H',p.Hs, 'JPD',p.JPD, ...
                                      'float_amplitude', val.X_f, ...
                                      'spar_amplitude',val.X_s, ...
                                      'relative_amplitude',val.X_u, ...
                                      'PTO_damping',val.B_p);
end

function results = load_wecsim_results(wecsim_filename, p)

    sz = size(p.JPD);
    % wecSim spar stationary
    vars = {'P','float_amplitude','spar_amplitude','relative_amplitude'};

    wecsim_raw = load(wecsim_filename,vars{:});
    power_mech_unsat = -reshape(wecsim_raw.P, sz) / 1000;
    %power(power>1e6) = NaN;
    float_amplitude = reshape(wecsim_raw.float_amplitude,sz);
    spar_amplitude = reshape(wecsim_raw.spar_amplitude,sz);
    relative_amplitude = reshape(wecsim_raw.relative_amplitude,sz);

    % todo: add damping


    results = assemble_results_struct(size(power_mech_unsat),...
                                         'power_mech_unsat',power_mech_unsat, ...
                                          'T', p.T, 'H', p.Hs, 'JPD', p.JPD, ...%PTO_damping
                                         'float_amplitude', float_amplitude,...
                                         'spar_amplitude',spar_amplitude, ...
                                         'relative_amplitude',relative_amplitude);

end

function results = assemble_results_struct(sz,varargin)

    var_names = {'power_mech_unsat', 'power_elec_unsat', 'power_elec_sat', ...
                 'T', 'H', 'JPD', 'float_amplitude', 'spar_amplitude', ...
                 'relative_amplitude', 'PTO_damping'};

    p = inputParser;
    for var_idx = 1:length(var_names)
        addParameter(p,var_names{var_idx},NaN(sz),@isnumeric);
    end
    parse(p,varargin{:})
    results = p.Results;

    % calculated variables
    wave_resource_raw = 1030 * 9.8^2 / (64*pi) * results.T .* results.H.^2 / 1000;
    %wave_resource_sim = wave_resource_raw;
    %wave_resource_sim(results.JPD == 0) = NaN;
    
    %diameter = X(2);
    results.CW = results.power_mech_unsat ./ wave_resource_raw;    % capture width
    CW_max = 9.8 * results.T.^2 / (4*pi^2);             % max CW for linear hydrodynamics, axisymmetric body
    %results.CWR = CW / diameter;                        % capture width ratio
    results.CW_to_CW_max = results.CW ./ CW_max;
    %CW_to_CW_max_zeroed = CW_to_CW_max;
    %JPD_rep = repmat(JPD,[1 1 4]);
    %CW_to_CW_max_zeroed(JPD_rep == 0) = NaN;

end

function pct_error = compute_percent_error_matrix(actual, sim)
    pct_error = 100 * (sim - actual) ./ actual;
    %pre_JPD_pct_error(JPD==0) = NaN;
    pct_error(pct_error==0) = NaN;

end