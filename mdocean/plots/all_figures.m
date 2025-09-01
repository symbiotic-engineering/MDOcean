% Generate all figures used in the paper
function [fig_success,tab_success,...
          fig_output, tab_output,...
          fig_runtime,tab_runtime,...
          num_figs,num_tabs,...
          fig_names,tab_names] = all_figures( which_figs, which_tabs, filename_uuid )

if nargin<3
    filename_uuid = ''; % required argument for anything running gradient_optim in parallel,
                        % since prob2struct needs unique filenames for code generation
end

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
save_folder = ['../test-results/' date '/'];
mkdir(save_folder)
table_save_fcn = @(tab,filename,varargin) table2latex(tab, [save_folder filename], varargin{:});

num_figs = 59;
num_tabs = 8;

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

fig_success = cell([1,length(which_figs)]);
tab_success = cell([1,length(which_tabs)]);

fig_output = gobjects(1, length(which_figs));
tab_output(1, 1:length(which_tabs)) = {table()};

fig_runtime = NaN([1,length(which_figs)]);
tab_runtime = NaN([1,length(which_tabs)]);
%%

figs_in_paper = cell([1,num_figs]);
figs_in_paper{1}  = 'ReadNonMatlabFigs.RM3_image';
% 2: modeling methodology
figs_in_paper{2}  = 'ReadNonMatlabFigs.methodology_overview';
figs_in_paper{3}  = 'ReadNonMatlabFigs.N2_diagram';
figs_in_paper{4}  = 'ReadNonMatlabFigs.dimensions';
figs_in_paper{5}  = 'ReadNonMatlabFigs.MEEM_geometry';
figs_in_paper{6}  = 'SparHydro.spar_added_mass';
figs_in_paper{7}  = 'HydroCoeffFigFunc.hydro_coeff_err';
figs_in_paper{8}  = 'DescFcns.saturation_desc_fcn';
figs_in_paper{9}  = 'DescFcns.drag_desc_fcn';
figs_in_paper{10} = 'RunSingleFigFunc.nominal_power_matrix';
figs_in_paper{11} = 'Wecsim.WECSim_error_histograms_multibody';
figs_in_paper{12} = 'ReadNonMatlabFigs.FBD';
figs_in_paper{13} = 'Validation.cost_vs_N_WEC';
figs_in_paper{14} = 'Runtime.sim_runtime';
figs_in_paper{15} = 'Runtime.hydro_runtime';
figs_in_paper{16} = 'Runtime.dynamics_runtime';
% 3: optimization methodology
figs_in_paper{17} = 'ReadNonMatlabFigs.optimization_flowchart';
% 4: results
figs_in_paper{18} = 'DesignSpaceExploration.experiments';
figs_in_paper{19} = 'Comparison.overlaid_geometry';
figs_in_paper{20} = 'Comparison.overlaid_hydro_coeffs';
figs_in_paper{21} = 'Comparison.probability_CDF';
figs_in_paper{22} = 'GradientOptimFigFunc.single_obj_opt_power_matrix';
figs_in_paper{23} = 'GradientOptimFigFunc.lagrange_multipliers';
figs_in_paper{24} = 'GradientOptimFigFunc.dJ_dx_gradient';
figs_in_paper{25} = 'Multistart.multistart_convergence_tree';
figs_in_paper{26} = 'Multistart.multistart_bar chart';
figs_in_paper{27} = 'ParamSensitivities.re_optim_objective_tornado';
figs_in_paper{28} = 'ParamSensitivities.re_optim_design_tornado_J1';
figs_in_paper{29} = 'ParamSensitivities.re_optim_design_tornado_J2';
figs_in_paper{30} = 'ParetoFigFunc.pareto_front_with_design_images';
figs_in_paper{31} = 'ParetoFigFunc.pareto_front_LCOE_contours';
figs_in_paper{32} = 'ParetoFigFunc.pareto_heuristics';
figs_in_paper{33} = 'ParetoFigFunc.pareto_constraint_activity';
figs_in_paper{34} = 'AllFigCompare.runtime_bar_chart';
figs_in_paper{35} = 'ParetoSweep.sweep_num_seeds';
% appendix A- hydro
figs_in_paper{36} = 'Meem.meem_regions';
figs_in_paper{37} = 'Meem.meem_sparsity';
figs_in_paper{38} = 'Meem.meem_validation';
figs_in_paper{39} = 'Meem.meem_matching';
figs_in_paper{40} = 'Meem.meem_convergence';
figs_in_paper{41} = 'Meem.asymptotic_b_vector';
% appendix B - dynamics
figs_in_paper{42} = 'ForceSaturationFigFunc.power_force_sensitivity';
figs_in_paper{43} = 'RunSingleFigFunc.drag_convergence';
figs_in_paper{44} = 'Slamming.slamming_amplitude';
figs_in_paper{45} = 'RunSingleFigFunc.slamming_model_comparison';
figs_in_paper{46} = 'Wecsim.wecsim_all_sea_states';
% appendix C - structures
figs_in_paper{47} = 'ReadNonMatlabFigs.equivalent_stiffness';
figs_in_paper{48} = 'ReadNonMatlabFigs.trapezoid';
figs_in_paper{49} = 'ReadNonMatlabFigs.damping_plate_flowchart';
figs_in_paper{50} = 'DampingPlateStructures.damping_plate_moment';
figs_in_paper{51} = 'DampingPlateStructures.damping_plate_deflection';
figs_in_paper{52} = 'DampingPlateStructures.damping_plate_aspect_ratio';
% appendix D - economics
% appendix E - optimization process
figs_in_paper{53} = 'ParamSensitivities.post_optim_re_optim_objective_grid';
figs_in_paper{54} = 'ParamSensitivities.post_optim_design_grid';
figs_in_paper{55} = 'ParamSensitivities.re_optim_design_grid';
% appendix F - supplementary results
figs_in_paper{56} = 'GradientOptimFigFunc.normalized_gradient';
figs_in_paper{57} = 'GradientOptimFigFunc.single_obj_convergence';
figs_in_paper{58} = 'Multistart.multistart_parallel_coordinates';
% graphical abstract (unnumbered so at the end)
figs_in_paper{59}  = 'ReadNonMatlabFigs.graphical_abstract';

%% TABLES
tabs_in_paper = cell(1,29);
tabs_in_paper{12} = 'Cost.cost_parameters';
tabs_in_paper{13} = 'Validation.validation';
tabs_in_paper{14} = 'Constraints.constraints';
tabs_in_paper{15} = 'DesignVars.design_vars';
tabs_in_paper{16} = 'Parameters.parameters';
tabs_in_paper{18} = 'Comparison.optimal_design_vars';
tabs_in_paper{19} = 'Comparison.optimal_outputs';
tabs_in_paper{20} = 'LocationSensitivity.location_sensitivity';

%% Auto-generate figure and table names from list above
fig_nums = cellstr(num2str((1:num_figs).'));
tab_idxs = ~cellfun(@isempty,tabs_in_paper);
tab_nums = cellstr(num2str(find(tab_idxs).'));

fig_names = strcat("Fig. ", fig_nums, ": ", figs_in_paper.');
tab_names = strcat("Tab. ", tab_nums, ": ", tabs_in_paper(tab_idxs).');

%% Initialize structures to store generated figures and tables
generated_figs = struct();
generated_tabs = struct();

% Create analysis class instances cache
analysis_instances = containers.Map();

% Loop for figures
for i = 1:length(which_figs)
    try
        % Get the class name for the figure
        fig_number = which_figs(i);
        fig_name = split(figs_in_paper{fig_number}, '.');
        class_name = fig_name{1};
        fig_desc = fig_name{2};
        
        % Check if the figure has already been generated
        t = tic;
        if ~isfield(generated_figs, class_name)
            % Create or retrieve analysis instance
            if ~isKey(analysis_instances, class_name)
                % Create new instance and set filename_uuid and table_save_fcn
                analysis_obj = eval([class_name '()']);
                analysis_obj.b.filename_uuid = filename_uuid;
                analysis_obj.b.table_save_fcn = table_save_fcn;
                analysis_instances(class_name) = analysis_obj;
            else
                analysis_obj = analysis_instances(class_name);
            end
            
            % Run the analysis if not already done
            if isempty(analysis_obj.fig_array)
                analysis_obj = analysis_obj.run_analysis();
                analysis_obj = analysis_obj.run_post_process();
                analysis_instances(class_name) = analysis_obj;
            end
            
            % Extract figures and tables from the analysis object
            figs = struct();
            tabs = struct();
            for j = 1:length(analysis_obj.fig_names)
                if j <= length(analysis_obj.fig_array)
                    figs.(analysis_obj.fig_names{j}) = analysis_obj.fig_array(j);
                end
            end
            for j = 1:length(analysis_obj.tab_names)
                if j <= length(analysis_obj.tab_array_display)
                    tabs.(analysis_obj.tab_names{j}) = analysis_obj.tab_array_display{j};
                end
            end
            
            % Store the generated figure in the generated_figs structure
            generated_figs.(class_name) = figs;
            
            % If the function generates a table, store it in the generated_tabs structure
            if ~isempty(tabs)
                generated_tabs.(class_name) = tabs;
            end
        end
        
        % Store the generated figure in fig_output
        fig_output(i) = generated_figs.(class_name).(fig_desc);
        fig_runtime(i) = toc(t);

    catch err
        fig_success{i} = err;  % Store error for the figure
    end
end

% Loop for tables
for i = 1:length(which_tabs)
    try
        tab_number = which_tabs(i);
        % Get the class name for the table
        tab_name = split(tabs_in_paper{tab_number}, '.');
        class_name = tab_name{1};
        tab_desc = tab_name{2};
        
        % Check if the table has already been generated
        t = tic;
        if ~isfield(generated_tabs, class_name)
            % Create or retrieve analysis instance
            if ~isKey(analysis_instances, class_name)
                % Create new instance and set filename_uuid and table_save_fcn
                eval(['analysis_obj = ' class_name '();']);
                analysis_obj.b.filename_uuid = filename_uuid;
                analysis_obj.b.table_save_fcn = table_save_fcn;
                analysis_instances(class_name) = analysis_obj;
            else
                analysis_obj = analysis_instances(class_name);
            end
            
            % Run the analysis if not already done
            if isempty(analysis_obj.tab_array_display)
                analysis_obj = analysis_obj.run_analysis();
                analysis_obj = analysis_obj.run_post_process();
                analysis_instances(class_name) = analysis_obj;
            end
            
            % Extract figures and tables from the analysis object
            figs = struct();
            tabs = struct();
            for j = 1:length(analysis_obj.fig_names)
                if j <= length(analysis_obj.fig_array)
                    figs.(analysis_obj.fig_names{j}) = analysis_obj.fig_array(j);
                end
            end
            for j = 1:length(analysis_obj.tab_names)
                if j <= length(analysis_obj.tab_array_display)
                    tabs.(analysis_obj.tab_names{j}) = analysis_obj.tab_array_display{j};
                end
            end
            
            % Store the generated table in the generated_tabs structure
            generated_tabs.(class_name) = tabs;
            
            % If the function generates a figure, store it in the generated_figs structure
            if ~isempty(figs)
                generated_figs.(class_name) = figs;
            end
        end
        
        % Store the generated table in tab_output
        tab_output{i} = generated_tabs.(class_name).(tab_desc);  % Store the table
        tab_runtime(i) = toc(t);
        display(tab_output{i})
        
    catch err
        tab_success{i} = err;  % Store error for the table
    end
end

end
