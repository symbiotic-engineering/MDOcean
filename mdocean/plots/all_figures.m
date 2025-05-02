% Generate all figures used in the paper
function [fig_success,tab_success,fig_output,tab_output,...
            num_figs,num_tabs,fig_names,tab_names] = all_figures( which_figs, which_tabs, filename_uuid )

if nargin<3
    filename_uuid = ''; % required argument for anything running gradient_optim in parallel,
                        % since prob2struct needs unique filenames for code generation
end

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
save_folder = ['../test-results/' date '/'];
mkdir(save_folder)

num_figs = 36;
num_tabs = 8;
fig_names   = cell([1,num_figs]);
tab_names   = cell([1,num_tabs]);

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

fig_success = cell([1,length(which_figs)]);
tab_success = cell([1,length(which_tabs)]);

fig_output = gobjects(1, length(which_figs));

tab_output(1, 1:length(which_tabs)) = {table()};
%%

figs_in_paper(1)  = 'read_non_matlab_figs.RM3_image';
% 2: modeling methodology
figs_in_paper(2)  = 'read_non_matlab_figs.methodology_overview';
figs_in_paper(3)  = 'read_non_matlab_figs.N2_diagram';
figs_in_paper(4)  = 'read_non_matlab_figs.dimensions';
figs_in_paper(5)  = 'read_non_matlab_figs.MEEM_geometry';
figs_in_paper(6) =  'spar_hydro.spar_added_mass';
figs_in_paper(7)  = 'hydro_coeff_err.hydro_coeff_err';
figs_in_paper(8)  = 'desc_fcns.saturation_desc_fcn';
figs_in_paper(9)  = 'desc_fcns.drag_desc_fcn';
figs_in_paper(10)  = 'run_single.power_matrix_multiply';
figs_in_paper(11) = 'wecsim.WECSim_error_histograms';
figs_in_paper(12) = 'read_non_matlab_figs.FBD';
figs_in_paper(13) = 'validation.cost_vs_N_WEC';
figs_in_paper(14) = 'runtime.sim_runtime';
figs_in_paper(15) = 'runtime.hydro_runtime';
figs_in_paper(16) = 'runtime.dynamics_runtime';
% 3: optimization methodology
figs_in_paper(17) = 'read_non_matlab_figs.optimization_flowchart';
% 4: results
figs_in_paper(18) = 'experiments.experiments';
figs_in_paper(19) = 'compare.overlaid_geometry';
figs_in_paper(20) = 'compare.overlaid_hydro_coeffs';
figs_in_paper(21) = 'compare.probability_CDF';
figs_in_paper(22) = 'gradient_optim.lagrange_multipliers';
figs_in_paper(23) = 'gradient_optim.dJ_dx_gradient';
figs_in_paper(24) = 'multistart.multistart_convergence_tree';
figs_in_paper(25) = 'param_sensitivities.objective_tornadoes';
figs_in_paper(26) = 'param_sensitivities.design_tornadoes';
figs_in_paper(27) = 'pareto.pareto_front_with_design_images';
figs_in_paper(28) = 'pareto.pareto_front_LCOE_contours';
figs_in_paper(29) = 'pareto.heuristics';
figs_in_paper(30) = 'pareto.constraint_activity';
figs_in_paper(31) = 'all_fig_compare.runtime_bar_chart';
% appendix A- hydro
figs_in_paper(32) = 'meem.meem_regions';
figs_in_paper(33) = 'meem.meem_sparsity';
figs_in_paper(34) = 'meem.meem_validation';
figs_in_paper(35) = 'meem.meem_matching';
figs_in_paper(36) = 'meem.meem_convergence';
figs_in_paper(37) = 'meem.asymptotic_b_vector';
% appendix B - dynamics
figs_in_paper(38) = 'run_single.drag_convergence';
figs_in_paper(39) = 'slamming.slamming_amplitude';
figs_in_paper(40) = 'force_sat_results.power_force_sensitivity';
figs_in_paper(41) = 'wecsim.wecsim_all_sea_states';
% appendix C - structures
figs_in_paper(42) = 'read_non_matlab_figs.trapezoid';
figs_in_paper(43) = 'read_non_matlab_figs.damping_plate_flowchart';
figs_in_paper(44) = 'damping_plate_structures.damping_plate_moment';
figs_in_paper(45) = 'damping_plate_structures.damping_plate_deflection';
figs_in_paper(46) = 'damping_plate_structures.damping_plate_aspect_ratio';
% appendix D - economics
% appendix E - optimization process
figs_in_paper(47) = 'param_sensitivities.local_global_objective';
figs_in_paper(48) = 'param_sensitivities.local_design';
figs_in_paper(49) = 'param_sensitivities.global_design';
% appendix F - supplementary results
figs_in_paper(50) = 'gradient_optim.single_obj_convergence';
% graphical abstract (unnumbered so at the end)
figs_in_paper(51)  = 'read_non_matlab_figs.graphical_abstract';

%% TABLES
tabs_in_paper = cell(1,29);
tabs_in_paper(12) = 'cost.cost_parameters';
tabs_in_paper(13) = 'validation.validation';
tabs_in_paper(14) = 'constraints';
tabs_in_paper(15) = 'design variables';
tabs_in_paper(16) = 'parameters';
tabs_in_paper(18) = 'compare.optimal_design';
tabs_in_paper(19) = 'compare.optimal_outputs';
tabs_in_paper(20) = 'location_sensitivity.location_sensitivity';

%%
fig_names{6} = 'Fig. 6: hydro coeffs vs freq';
fig_names{7} = 'Fig. 7: drag describing function';
fig_names{8} = 'Fig. 8: force saturation time signal';
fig_names{9} = 'Fig. 9: force saturation results';
fig_names{10} = 'Fig. 10: JPD multiplication';
fig_names{12} = 'Fig. 12: cost vs N WEC';
fig_names{15} = 'Fig. 15: design space exploration';
fig_names{16} = 'Fig. 16: local and global objective parameter sensitivities';
fig_names{17} = 'Fig. 17: local optimal design variable parameter sensitivities';
fig_names{18} = 'Fig. 18: global optimal design variable parameter sensitivities';
fig_names{22} = 'Fig. 22: pareto front with design images';
fig_names{23} = 'Fig. 23: design and objective heuristics';
fig_names{24} = 'Fig. 24: pareto front with LCOE contours';
fig_names{25} = 'Fig. 25: constraint activity';
fig_names{30} = 'Fig. 30: asymptotic b vector';

tab_names{1} = 'Tab. 12: validation against nominal';
tab_names{2} = 'Tab. 15: constraints';
tab_names{3} = 'Tab. 16: design variables';
tab_names{4} = 'Tab. 17: parameters';
tab_names{5} = 'Tab. 19: optimal DVs for various designs';
tab_names{6} = 'Tab. 20: optimal outputs for various designs';
tab_names{7} = 'Tab. 21: convergence for different x0s';
tab_names{8} = 'Tab. 22: optimal DVs for 4 locations';

%% Initialize structures to store generated figures and tables
generated_figs = struct();
generated_tabs = struct();
p = parameters();
b = var_bounds();
b.filename_uuid = filename_uuid;

% Loop for figures
for i = which_figs
    try
        % Get the function name for the figure

        fig_name = split(figs_in_paper{i}, '.');
        func_name = fig_name{1};
        fig_desc = fig_name{2};
        func_handle = str2func(func_name);
        
        % Check if the figure has already been generated
        if ~isfield(generated_figs, func_name)
            % Run the function to generate the figure and table
            [figs, tabs] = func_handle(p, b);
            
            % Store the generated figure in the generated_figs structure
            generated_figs.(func_name) = figs;
            
            % If the function generates a table, store it in the generated_tabs structure
            if ~isempty(tabs)
                generated_tabs.(func_name) = tabs;
            end
        end
        
        % Store the generated figure in fig_output
        fig_output(i) = generated_figs.(func_name).(fig_desc);
        
    catch err
        fig_success{i} = err;  % Store error for the figure
    end
end

% Loop for tables
for i = which_tabs
    try
        % Get the function name for the table
        tab_name = split(tabs_in_paper{i}, '.');
        func_name = tab_name{1};
        tab_desc = tab_name{2};
        func_handle = str2func(func_name);
        
        % Check if the table has already been generated
        if ~isfield(generated_tabs, func_name)
            % Run the function to generate the table
            [figs, tabs] = func_handle(p, b);
            
            % Store the generated table in the generated_tabs structure
            generated_tabs.(func_name) = tabs;
            
            % If the function generates a figure, store it in the generated_figs structure
            if ~isempty(figs)
                generated_figs.(func_name) = figs;
            end
        end
        
        % Store the generated table in tab_output
        tab_output{i} = generated_tabs.(func_name).(tab_desc);  % Store the table
        display(tab_output{i})
        
    catch err
        tab_success{i} = err;  % Store error for the table
    end
end

end

%%
% Non-MATLAB figure reading function
function [figs,tabs] =  read_non_matlab_figs_func(~,~)

    fig_names = {'RM3_image','methodology_overview','N2_diagram',...
                'dimensions','MEEM_geometry','WECSim_error_breakdown', ...
                'FBD', 'sim_runtime', 'optimization_flowchart'};
    file_names = {'geometry.png',...
                    'methods_flowchart_2_cropped.jpg',...
                    'xdsm.JPG',...
                    'dimensions.jpg',...
                    'MEEM-dims-basic-2.jpg',...
                    'Error_Accumulation_AEP.png',...
                    'structures_FBDs_4.jpg',...
                    'time_breakdown_3.png',...
                    'optim_process_flowchart_cropped.jpg'};
    % Initialize an empty struct to hold the figures
    figs = struct();
    tabs = [];

    % Loop over the requested figure names
    for i = 1:length(fig_names)
        fig_name = fig_names{i};
        file_name = file_names{i};
        
        % Create the figure and read in the image file
        figs.(fig_name) = figure;
        imshow(imread(file_name), 'Parent', axes(figs.(fig_name)));
        figs.(fig_name).UserData = ['plots/non_matlab_figs/' file_name(1:end-4) '.pdf'];
    end
end

% Function for single objective convergence
function [figs,tabs] = gradient_optim_func(p,b)
    gradient_optim(b.X_start_struct,p,b,1,{@optimplotfval, @(x,~,~)optim_geomviz(x,p,b)})
    figs.convergence = gcf;
end

% Function to generate Pareto figure group
function [figs,tabs] =  pareto_fig_func(~,~)
    pareto_search(filename_uuid)
    pareto_design_heuristics()
    figs.pareto_heuristics = gcf;
    n = figs.pareto_heuristics.Number;
    figs.pareto_with_design_images = figure(n - 2);
    figs.pareto_front_LCOE_contours = figure(n - 3);
    figs.pareto_constraint_activity = figure(n - 5); 
    figs.pareto_constraint_activity.Position = [1 41 1536 844.8000];
    tabs = [];
end

% Function to generate parameter sensitivity figures
function [figs,tabs] =  param_sensitivities_fig_func(~,b)
    param_sweep(filename_uuid)
    figs.global_constraint = gcf;                   % delta p global (grid)
    n = figs.global_constraint.Number;
    figs.local_constraint= figure(n-1);             % delta p local (grid)
    figs.global_design = figure(n - 2);             % dx*/dp global (grid)
    figs.local_design  = figure(n - 3);             % dx*/dp local (grid)
    figs.global_objective = figure(n - 4);          % dJ*/dp global (grid)
    figs.local_global_objective = figure(n - 5);    % dJ*/dp combined (grid)
    figs.global_design_tornado_J2 = figure(n - 6);  % dx*/dp global (tornado)
    figs.global_design_tornado_J1 = figure(n - 7);  % dx*/dp global (tornado)
    num_DVs = length(b.var_names);
    for i = 1:num_DVs
        name = ['global_design_tornado_' b.var_names(num_DVs+1-i)];
        figs.(name) = figure(n - (7+i));            % dx*/dp global (tornado)
    end
    figs.nonlinear_design_J2 = figure(n - num_DVs - 9);
    figs.nonlinear_design_J1 = figure(n - num_DVs - 10);
    figs.nonlinear_objectives = figure(n - num_DVs - 11);
    tabs = [];
end

% Function to generate hydro coefficient figure
function [figs,tabs] =  hydro_coeff_fig_func(~,~)
    hydro_coeff_err()
    figs.hydro_coeff_err = gcf;
    tabs = [];
end

% Function to generate drag describing function and saturation time signal figures
function [figs,tabs] =  drag_saturation_fig_func(~,~)
    sin_desc_fcn_demo()
    figs.drag_desc_fcn = gcf;
    figs.saturation_desc_fcn = figure(figs.drag_desc_fcn.Number - 2);
    tabs = [];
end

% Function to generate force saturation results figure
function [figs,tabs] =  force_saturation_fig_func(p, b)
    fig_array = force_sat_results(p, b);
    figs.power_force_sensitivity = fig_array(1);
    figs.runtime_sensitivity = fig_array(2);
    tabs = [];
end

% Function to generate JPD multiplication figure
function [figs,tabs] = jpd_multiply_fig_func(p, b)
    X = [b.X_noms; 1];
    plot_power_matrix(X, p, filename_uuid)
    figs.JPD_multiplication = gcf;
    tabs = [];
end

% Function to generate validation figure for cost vs N WEC
function [figs,tabs] = cost_vs_nwec_fig_func(~,~)

    [~, ~, ~, ~, tab1a, figs.cost_vs_N_WEC] = validate_nominal_RM3('report');
    [~,~,~,~,tab1b] = validate_nominal_RM3('wecsim');

    % merge table 1a and 1b while preserving row order
    sharedCols = intersect(tab1a.Properties.VariableNames, tab1b.Properties.VariableNames);
    tab1a.RowNum = (1:length(tab1a.Properties.RowNames))';
    tab1b.RowNum = length(tab1a.Properties.RowNames) + (1:length(tab1b.Properties.RowNames))';
    
    tab1 = outerjoin(tab1a, tab1b, 'Keys', [{'RowNum'},sharedCols], 'MergeKeys', true);
    tab1 = removevars(tab1,'RowNum');
    tab1.Properties.RowNames = [tab1a.Properties.RowNames; tab1b.Properties.RowNames];

    vector_cols = {'capex','capex_design','J_capex_design','capex_struct','capex_PTO','opex','LCOE'};
    idx_remove = ismember(tab1.Properties.VariableNames,vector_cols); 
    tab1latex = rows2vars(tab1,'VariableNamingRule','preserve','DataVariables',~idx_remove);
    tab1latex.OriginalVariableNames = remove_underscores(modify_suffix(tab1latex.OriginalVariableNames));
    new_names = {'Variable','MDOcean','Actual','Error','MDOcean ','Actual ','Error '};
    tab1latex = renamevars(tab1latex, tab1latex.Properties.VariableNames, new_names);
    firstrow = '&\multicolumn{3}{c|}{DOE Report RM3 Design \cite{RM3}} & \multicolumn{3}{c}{WEC-Sim RM3 Design} \\';
    colspec = '>{\centering\arraybackslash}p{0.2\linewidth}|c|c|r|c|c|r';
    table2latex(tab1latex,[save_folder 'table_12.tex'],colspec,firstrow)

    tabs.validation_table = tab1;
end

% Function to generate design space exploration figure
function [figs,tabs] = design_space_exploration_fig_func(p,b)
    experiments(p,b)
    figs.experiments = gcf;
    tabs = [];
end

% Function to generate overlaid geometry, hydro coeffs, and probability CDF figures
function [figs,tabs] = comparison_fig_tab_func(p,b)
    [tabs.optimal_design_vars, tabs.optimal_outputs] = compare(p,b);
    n = gcf().Number;
    figs.overlaid_geometry = figure(n-3);
    figs.overlaid_hydro_coeffs = figure(n);
    figs.probability_CDF = figure(n-2);

    table2latex(tabs.optimal_design_vars,[save_folder 'table_19.tex'])
    table2latex(tabs.optimal_outputs,    [save_folder 'table_20.tex'])
end

% Function to generate asymptotic b vector figure
function [figs,tabs] = asymptotic_b_vector_fig_func(~,~)
    b_inf_numeric()
    figs.asymptotic_b_vector = gcf;
    tabs = [];
end

% Function to generate damping plate moment, deflection, and aspect ratio figures
function [figs,tabs] = damping_plate_fig_func(~,~)
    BoedoPrantilAnnularPlate()
    n = gcf().Number;
    figs.damping_plate_aspect_ratio = figure(n);
    figs.damping_plate_deflection   = figure(n - 1);
    figs.damping_plate_moment       = figure(n - 7);
    tabs = [];
end

function [figs,tabs] = location_sensitivity_func(p,b)
    tab = location_sensitivity(p,b);
    tabs.location_table = tab;
    location_flags = str2double(tab(strcmp(tab.Row,'flag'),:).Variables);
    
    %tab_success{which_tabs == 8} = {location_flags};

    idx_remove = ismember(tab.Row,{'flag','Optimal Material index'});
    tablatex = tab(~idx_remove,:);
    colspec = ['>{\centering\arraybackslash}p{0.20\linewidth}' ...
               '>{\centering\arraybackslash}p{0.08\linewidth}' ...
               '>{\centering\arraybackslash}p{0.15\linewidth}' ...
               '>{\centering\arraybackslash}p{0.15\linewidth}' ...
               '>{\centering\arraybackslash}p{0.15\linewidth}' ...
               '>{\centering\arraybackslash}p{0.18\linewidth}'];
    firstrow = '&& \multicolumn{4}{c}{Location}\\  \cline{3-6}';
    table2latex(tablatex,[save_folder 'table_22.tex'],colspec,firstrow)

    figs = [];
end

function [figs,tabs] = constraints_func(~,b)
    tab = b.constraint_names';
    tabs.constraints = tab;
    figs = [];
end

function [figs,tabs] = design_vars_func(~,b)
    tab = array2table([b.X_mins b.X_noms b.X_maxs], ...
            'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', b.var_names(1:end-1));
    tabs.design_vars = tab;
    figs = [];
end

function [figs,tabs] = params_func(~,~)
    [~,tabs.parameters] = parameters();
    figs = [];
end

function [figs,tabs] = multistart(p,b)
    tab = gradient_mult_x0(p,b);
    tabs.multistart = tab;
    table2latex(tab,[save_folder 'table_21.tex'])
    figs = [];
end

% %% figure 14-15 - wecsim validation histograms
% fig_names{14} = 'Fig. 14: WecSim histogram singlebody';
% fig_names{15} = 'Fig. 15: WecSim histogram multibody';
% if any(which_figs == 14 | which_figs == 15)
%     try
%         [~, ~, ~, ~, fig_singlebody, fig_multibody] = validate_dynamics();
%         fig_output(which_figs==14) = fig_singlebody;
%         fig_output(which_figs==15) = fig_multibody;
%     catch err
%         fig_success(which_figs == 14 | which_figs == 15) = {err};
%     end
% end


%% figure 31-33 - fixme not implemented
if any( ismember(which_figs, 31:33) )
    fig_output(ismember(which_figs, 31:33)) = figure;
end

%% fixme: 37 to 45 not implemented, but still name the figures
empty_idx = cellfun(@isempty,fig_names);
empty_str = strcat('Fig._', string(find(empty_idx)));
[fig_names{empty_idx}] = deal(empty_str{:});
[fig_success{ismember(which_figs,find(empty_idx))}] = deal(MException('all_figures:not_implemented','Figure not implemented yet'));


end
