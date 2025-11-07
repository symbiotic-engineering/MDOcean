clear
close all

file_95_reactive = 'C:\Users\rgm222\Downloads\test-report(95)\mdocean\wecsim_sparcd0_floatcd0_multibody_0_meem_off_git_746f2f5_uuid_ed37317d-3b76-4d02-9ea4-480325353cf8';
file_96_damping = 'C:\Users\rgm222\Downloads\test-report(96)\mdocean\wecsim_sparcd0_floatcd0_multibody_0_meem_off_git_ef9a674_uuid_4f19928e-8e18-434d-b7c7-96d67ba98c60';

reac = load(file_95_reactive);
damp = load(file_96_damping);

results_actual = load_wecsim_results(file_95_reactive, reac.p);
results_sim    = load_wecsim_results(file_96_damping,  damp.p);

vars_to_plot = {'power_mech_unsat','power_elec_sat','CW_to_CW_max',...
        'float_amplitude','relative_amplitude','spar_amplitude','PTO_damping','force_pto'};
figs = comparison_plot(reac.p.T, reac.p.Hs, results_actual, results_sim, vars_to_plot, {'reactive'}, {'damping'}, reac.p);


%% copy pasted from power matrix compare
function weighted_error = compute_weighted_percent_error(sim, actual, weights)
    sim_weighted = sim .* weights / 100;
    sim_avg = sum(sim_weighted(:),'omitnan');

    actual_weighted = actual .* weights / 100;
    actual_avg = sum(actual_weighted(:),'omitnan');
    weighted_error = compute_percent_error_matrix(actual_avg, sim_avg);
end

function figs = comparison_plot(T, H, actual, sim, vars_to_plot, actual_str, sim_str, p)
    figs = gobjects(1,length(vars_to_plot));
    num_subplots = 2*length(sim)+1;
    for var_idx = 1:length(vars_to_plot)
        var_name = vars_to_plot{var_idx};

        % find max and min to use as consistent colorbar axis for sim and actual
        flatten = @(a) a(:);
        all_data = [flatten(actual.(var_name)); flatten([sim(:).(var_name)])];
        clims = [min(all_data,[],'all'),max(all_data,[],'all')];
        if all(clims==[0 0]) || all(isnan(clims))
            clims = [0 1];
        end

        f = figure;
        figs(var_idx) = f;
        % actual
        subplot(1,num_subplots,1)
        contour_plot(T,H, actual.(var_name), ['Actual: ',actual_str]);
        caxis(clims)

        for i=1:length(sim)
            % sim
            subplot(1,num_subplots,1+i)
            contour_plot(T,H,sim(i).(var_name), ['Sim: ',sim_str{i}] );
            caxis(clims)
            % error
            error = compute_percent_error_matrix(actual.(var_name), sim(i).(var_name));
            subplot(1,num_subplots,1+length(sim)+i)

            % set contour lines to make plot more readable
            error_levels = [];
            if strcmp(var_name,'power_mech_unsat') || strcmp(var_name,'CW_to_CW_max')
                if     p.C_d_float==0 && p.use_MEEM==false && p.use_multibody==false
                    error_levels = 30:2:48;
                elseif p.C_d_float==1 && p.use_MEEM==false && p.use_multibody==false
                    error_levels = [-10 0 10 35 50 100 150 200 250 300 350 400];
                elseif p.C_d_float==0 && p.use_MEEM==true  && p.use_multibody==false
                    error_levels = [-25:5:0 2 5];
                elseif p.C_d_float==1 && p.use_MEEM==true  && p.use_multibody==false
                    error_levels = [-85 -50 -20 -14 -12 -10:5:20];
                end
            elseif strcmp(var_name,'float_amplitude') || strcmp(var_name,'relative_amplitude')
                if     p.C_d_float==0 && p.use_MEEM==false && p.use_multibody==false
                    error_levels = [-5:1:0 1.5 3];
                elseif p.C_d_float==1 && p.use_MEEM==false && p.use_multibody==false
                    error_levels = [-50 -40 -20 -10 0:5 10];
                elseif p.C_d_float==0 && p.use_MEEM==true  && p.use_multibody==false
                    error_levels = -20:5:10;
                elseif p.C_d_float==1 && p.use_MEEM==true  && p.use_multibody==false
                    error_levels = [-100 -10 -8 -7 -5 -2 0 2 5 10 20 30];
                elseif p.C_d_float==0 && p.use_MEEM==true && p.use_multibody==true
                    error_levels = -10:5:40;
                end
            end

            error_plot(T,H,error,['Percent Error ' sim_str{i}],error_levels);
        end

        sgtitle(remove_underscores(var_name))
    end
end

function error_plot(T, H, error, error_title, error_values)
    [c,h_fig] = contour_plot(T, H, error, error_title, error_values);
    clabel(c,h_fig);
    if ~isempty(h_fig)
        colormap(h_fig.Parent,bluewhitered)
    end
end

function [c,h_fig] = contour_plot(T, H, Z, Z_title, Z_levels)
    if nargin<5
        Z_levels = [];
    end

    if any(isfinite(Z),'all')
        [c,h_fig] = contourf(T,H,Z);
    else
        c = []; h_fig = [];
    end
    title(Z_title)
    xlabel('Wave Period T (s)')
    ylabel('Wave Height Hs (m)')
    cb = colorbar;
    if ~isempty(h_fig)
        levs = sort([cb.Ticks min(Z,[],'all') max(Z,[],'all') Z_levels]);
        h_fig.LevelList = levs;
    end
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

    JPD_temp = readmatrix(report_filename,'Range','E24:S37','Sheet',sheet);
    JPD = JPD_temp/sum(JPD_temp,'all');

    wave_resource_sheet = readmatrix(report_filename,'Range','E49:S62','Sheet',sheet);
    wave_resource_sheet(wave_resource_sheet == 0) = NaN;

    if override
        power_mech_unsat = power_mech_unsat(3:4,4:5);
        power_elec_sat = power_elec_sat(3:4,4:5);
        Hs = Hs(3:4); Te = Te(4:5);
        JPD = JPD(3:4,4:5);
        wave_resource_sheet = wave_resource_sheet(3:4,4:5);
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
    X_unsat = X;
    idx_P = 7;
    X_unsat(idx_P) = 1e9;
    [~, P_matrix, ~, val] = simulation(X_unsat,p);
    power_elec_unsat = P_matrix/1000;
    power_mech_unsat = val.P_mech/1000;

    % saturated power
    [~, P_matrix] = simulation(X,p);
    power_elec_sat = P_matrix/1000;

    results = assemble_results_struct(size(power_mech_unsat),...
                                      'power_mech_unsat', power_mech_unsat, ...
                                      'power_elec_unsat',power_elec_unsat, ...
                                      'power_elec_sat', power_elec_sat, ...
                                      'T',p.T, 'H',p.Hs, 'JPD',p.JPD, ...
                                      'float_amplitude', val.X_f, ...
                                      'spar_amplitude',val.X_s, ...
                                      'relative_amplitude',val.X_u, ...
                                      'PTO_damping',val.B_p,...
                                      'force_pto',val.mag_U);
end

function results = load_wecsim_results(wecsim_filename, p)

    sz = size(p.JPD);
    % wecSim spar stationary
    vars = {'P','float_amplitude','spar_amplitude','relative_amplitude','force_pto','B_p'};
    file = 'Humboldt_California_Wave Resource _SAM CSV.csv';
    old_jpd = trim_jpd(readmatrix(file,'Range','A3'));
    JPD = p.JPD; % old_jpd(2:end,2:end) - uncomment to use old wec sim results locally
    idx = find(JPD ~= 0);

    wecsim_raw = load(wecsim_filename,vars{:});
    [power_mech_unsat,...
     float_amplitude,...
     spar_amplitude,...
     relative_amplitude,...
     force_pto,...
     B_p] = deal(nan(sz));

    power_mech_unsat(idx) = -wecsim_raw.P / 1000;
    float_amplitude(idx) = wecsim_raw.float_amplitude;
    spar_amplitude(idx) = wecsim_raw.spar_amplitude;
    relative_amplitude(idx) = wecsim_raw.relative_amplitude;
    force_pto(idx) = wecsim_raw.force_pto;
    B_p(idx) = wecsim_raw.B_p;

    %power_mech_unsat(power_mech_unsat>1e6) = NaN;

    results = assemble_results_struct(size(power_mech_unsat),...
                                         'power_mech_unsat',power_mech_unsat, ...
                                          'T', p.T, 'H', p.Hs, 'JPD', p.JPD, ...
                                          'PTO_damping', B_p,'force_pto',force_pto,...
                                         'float_amplitude', float_amplitude,...
                                         'spar_amplitude',spar_amplitude, ...
                                         'relative_amplitude',relative_amplitude);

end

function results = assemble_results_struct(sz,varargin)

    var_names = {'power_mech_unsat', 'power_elec_unsat', 'power_elec_sat', ...
                 'T', 'H', 'JPD', 'float_amplitude', 'spar_amplitude', ...
                 'relative_amplitude', 'PTO_damping','force_pto'};

    p = inputParser;
    for var_idx = 1:length(var_names)
        addParameter(p,var_names{var_idx},NaN(sz),@isnumeric);
    end
    parse(p,varargin{:})
    results = p.Results;

    % calculated variables
    wave_resource_raw = 1030 * 9.8^2 / (64*pi) * results.T .* results.H.^2 / 1000;

    %diameter = X(2);
    results.CW = results.power_mech_unsat ./ wave_resource_raw;    % capture width
    CW_max = 9.8 * results.T.^2 / (4*pi^2);             % max CW for linear hydrodynamics, axisymmetric body
    %results.CWR = CW / diameter;                        % capture width ratio
    results.CW_to_CW_max = results.CW ./ CW_max;

end

function pct_error = compute_percent_error_matrix(actual, sim)
    pct_error = 100 * (sim - actual) ./ actual;
    pct_error(pct_error==0) = NaN;

end
