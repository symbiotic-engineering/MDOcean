% Generate all figures used in the paper
function [success_criterion,fig_output,tab_output,...
            num_figs,num_tabs,fig_names,tab_names] = all_figures( which_figs, which_tabs, filename_uuid )

if nargin<3
    filename_uuid = ''; % required argument for anything running gradient_optim in parallel,
                        % since prob2struct needs unique filenames for code generation
end

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
save_folder = ['../test-results/' date '/'];
mkdir(save_folder)

num_figs = 37;
num_tabs = 8;
fig_names = cell([1,num_figs]);
tab_names = cell([1,num_tabs]);
success_criterion = {};

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

if isempty(which_figs)
    which_figs = 0;
end
if isempty(which_tabs)
    which_tabs = 0;
end

fig_output = gobjects(1, length(which_figs));

tab_output(1, 1:length(which_tabs)) = {table()};

%% Define mapping from figures/tables to scripts
% figure_mapping = zeros(1,num_figs);
% figure_mapping(3:4) = 1; % sin desc function demo
% figure_mapping(5) = 2; % JPD
% figure_mapping(6:7) = 3; % pareto front
% figure_mapping(8) = 4; % parameter sensitivities
% figure_mapping(9:10) = 5; % overlaid comparison
% 
% table_mapping = zeros(1,num_tabs);
% table_mapping(1:2) = 6; % var bounds
% table_mapping(3) = 7; % parameters
% table_mapping(4) = 8; % validation
% table_mapping(5) = 5; % 3 design comparison
% table_mapping(6) = 9; % location sensitivity
% table_mapping(7) = 10; % convergence


%% non matlab figures
non_matlab_figs = [1:5 10:11 13:14];
non_matlab_fig_names = {'Fig. 1: RM3 image',...
                        'Fig. 2: Methodology overview',...
                        'Fig. 3: N2 diagram',...
                        'Fig. 4: Dimensions',...
                        'Fig. 5: MEEM geometry',...
                        'Fig. 10: WECSim error breakdown',...
                        'Fig. 11: FBD',...
                        'Fig. 13: sim runtime',...
                        'Fig. 14: Optimization flowchart'};
non_matlab_fig_files = {'geometry.png',...
                        'methods_flowchart_2_cropped.jpg',...
                        'xdsm.JPG',...
                        'dimensions.jpg',...
                        'MEEM-dims-basic-2.jpg',...
                        'Error_Accumulation_AEP.png',...
                        'structures_FBDs_4.jpg',...
                        'time_breakdown_3.png',...
                        'optim_process_flowchart_cropped.jpg'};

for i = 1:length(non_matlab_figs)
    fig_num = non_matlab_figs(i);
    fig_name = non_matlab_fig_names(i);
    fig_file = non_matlab_fig_files{i};

    fig_names{fig_num} = fig_name;
    if any(which_figs == fig_num)
        tmp_fig = figure;
        imshow(imread(fig_file),'Parent',axes(tmp_fig));
        tmp_fig.UserData = ['plots/non_matlab_figs/' fig_file(1:end-4) '.pdf'];
        fig_output(which_figs==fig_num) = tmp_fig;
    end
end

%% figure 6 - hydro coeffs vs freq
fig_names{6} = 'Fig. 6: hydro coeffs vs freq';
if any(which_figs == 6)
    hydro_coeff_err()
    fig6 = gcf;
    fig_output(which_figs==6) = fig6;
end

%% figure 7, 8 - drag DF, saturation time signal
fig_names{7} = 'Fig. 7: drag describing function';
fig_names{8} = 'Fig. 8: force saturation time signal';
if any(which_figs == 7 | which_figs == 8)
    sin_desc_fcn_demo()
    fig7 = gcf;
    fig8 = figure(fig7.Number-2);
    fig_output(which_figs==7) = fig7;
    fig_output(which_figs==8) = fig8;
end

%% figure 9 - JPD multiplication
fig_names{9} = 'Fig. 9: JPD multiplication';
if any(which_figs == 9)
    p = parameters();
    b = var_bounds();
    X = [b.X_noms; 1];
    plot_power_matrix(X,p)
    fig9 = gcf;
    fig_output(which_figs==9) = fig9;
end

%% figure 12 - cost vs N WEC
fig_names{12} = 'Fig. 12: cost vs N WEC';
if any(which_figs == 12) || any(which_tabs == 1)
    [~,~,~,~,tab1a,fig12] = validate_nominal_RM3('report');
    fig_output(which_figs==12) = fig12;
end

%% figure 15 - design space exploration
fig_names{15} = 'Fig. 15: design space exploration';
if any(which_figs == 15)
    % fixme - not implemented
    fig15 = figure;
    fig_output(which_figs==15) = fig15;
end

%% figure 16, 17, 18, 19, 20, 21 - parameter sensitivities
fig_names{16} = 'Fig. 16: local objective parameter sensitivities';
fig_names{17} = 'Fig. 17: local optimal objective parameter sensitivities';
fig_names{18} = 'Fig. 18: local optimal design variable parameter sensitivities';
fig_names{19} = 'Fig. 19: constraint activity thresholds';
fig_names{20} = 'Fig. 20: global optimal objective parameter sensitivities';
fig_names{21} = 'Fig. 21: global optimal design variable parameter sensitivities';
if any(which_figs == 16)
    param_sweep(filename_uuid)
    fig16 = gcf;
    fig_output(which_figs==16) = fig16;
    % fixme: 17 through 21 not implemented
end

%% figure 22, 23, 24, 25, 26 - pareto front, design heuristics
fig_names{22} = 'Fig. 22: pareto front with design images';
fig_names{23} = 'Fig. 23: design heuristics';
fig_names{24} = 'Fig. 24: objective heuristics';
fig_names{25} = 'Fig. 25: pareto front with LCOE contours';
fig_names{26} = 'Fig. 26: constraint activity';
if any(which_figs == 22 | which_figs == 23 | which_figs == 24 | which_figs == 25 | which_figs == 26)  
    pareto_search(filename_uuid);
    pareto_curve_heuristics()
    fig24 = gcf;
    fig23 = figure(fig24.Number - 1);
    fig22 = figure(fig24.Number - 3);
    fig25 = figure(fig24.Number - 4);
    fig26 = figure(fig24.Number - 6); % constraint activity
    fig26.Position = [1 41 1536 844.8000];
    fig_output(which_figs==22) = fig22;
    fig_output(which_figs==23) = fig23;
    fig_output(which_figs==24) = fig24;
    fig_output(which_figs==25) = fig25;
    fig_output(which_figs==26) = fig26;
end

%% figure 27, 28, 29 - overlaid geometry, hydro coeffs, probability CDF
fig_names{27} = 'Fig. 27: overlaid geometry';
fig_names{28} = 'Fig. 28: overlaid hydro coeffs';
fig_names{29} = 'Fig. 29: probability CDF';
if any(which_figs == 27 | which_figs == 28 | which_figs == 29 | which_tabs == 5 | which_tabs == 6)
    [tab5,tab6] = compare(filename_uuid);
    n = gcf().Number;
    fig29 = figure(n-1);
    fig27 = figure(n-2);
    fig_output(which_figs==27) = fig27;
    fig_output(which_figs==28) = figure;
    fig_output(which_figs==29) = fig29;
    % fixme: 28 not implemented
end

%% figure 30 - fixme not implemented
fig_names{30} = 'Fig. 30: XXX';
if any(which_figs == 30)
    fig_output(which_figs==30) = figure;
end

%% figure 31 - asymptotic b vector
fig_names{31} = 'Fig. 31: asymptotic b vector';
if any(which_figs == 31)
    b_Inf
    fig31 = gcf;
    fig_output(which_figs==31) = fig31;
end

%% figure 32-36 - fixme not implemented
if any( ismember(which_figs, 32:36) )
    fig_output(ismember(which_figs, 32:36)) = figure;
end

%% figure 37 - damping plate moment
fig_names{37} = 'Fig. 37: damping plate moment';
if any(which_figs == 37)
    BoedoPrantilAnnularPlate()
    fig37 = gcf;
    fig_output(which_figs==37) = fig37;
end

% fixme: 38 to 46 not implemented, but still name the figures
empty_idx = cellfun(@isempty,fig_names);
empty_str = strcat('Fig._', string(find(empty_idx)));
[fig_names{empty_idx}] = deal(empty_str{:});

%% table 12 - validation table
tab_names{1} = 'Tab. 12: validation against nominal';
if any(which_tabs == 1)
    % tab1a was generated with fig 12 above
    display(tab1a)
    [~,~,~,~,tab1b] = validate_nominal_RM3('wecsim');
    display(tab1b)

    % merge table 1a and 1b while preserving row order
    sharedCols = intersect(tab1a.Properties.VariableNames, tab1b.Properties.VariableNames);
    tab1a.RowNum = (1:length(tab1a.Properties.RowNames))';
    tab1b.RowNum = length(tab1a.Properties.RowNames) + (1:length(tab1b.Properties.RowNames))';
    
    tab1 = outerjoin(tab1a, tab1b, 'Keys', [{'RowNum'},sharedCols], 'MergeKeys', true);
    tab1 = removevars(tab1,'RowNum');
    tab1.Properties.RowNames = [tab1a.Properties.RowNames; tab1b.Properties.RowNames];
    display(tab1)

    tab_output{which_tabs==1} = tab1;
end

%% table 15 - constraints table
tab_names{2} = 'Tab. 15: constraints';
if any(which_tabs == 2)
    b = var_bounds();
    tab2 = b.constraint_names';
    display(tab2)
    tab_output{which_tabs==2} = tab2;
end

%% table 16 - design variables table
tab_names{3} = 'Tab. 16: design variables';
if any(which_tabs == 3)
    b = var_bounds();
    tab3 = array2table([b.X_mins b.X_noms b.X_maxs], ...
        'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', b.var_names(1:end-1));
    display(tab3)
    tab_output{which_tabs==3} = tab3;
end

%% table 17 - parameters table
tab_names{4} = 'Tab. 17: parameters';
if any(which_tabs == 4)
    [~,tab4] = parameters();
    display(tab4)
    tab_output{which_tabs==4} = tab4;
end

%% table 19 - optimal DVs for various designs
tab_names{5} = 'Tab. 19: optimal DVs for various designs';
if any(which_tabs == 5)
    % computation above with figures 27-29
    display(tab5);
    tab_output{which_tabs==5} = tab5;
    table2latex(tab5,[save_folder 'table_19.tex'])
end

%% table 20 - optimal outputs for various designs
tab_names{6} = 'Tab. 20: optimal outputs for various designs';
if any(which_tabs == 6)
    % computation above with figures 27-29
    display(tab6);
    tab_output{which_tabs==6} = tab6;
    table2latex(tab6,[save_folder 'table_20.tex'])
end

%% table 21 - convergence for different x0s
tab_names{7} = 'Tab. 21: convergence for different x0s';
if any(which_tabs == 7)
    tab7 = gradient_mult_x0(filename_uuid);
    tab_output{which_tabs==7} = tab7;
    table2latex(tab7,[save_folder 'table_21.tex'])
end

%% table 22 - optimal DVs for 4 locations
tab_names{8} = 'Tab. 22: optimal DVs for 4 locations';
if any(which_tabs == 8)
    tab8 = location_sensitivity(filename_uuid);
    display(tab8);
    tab_output{which_tabs==8} = tab8;
    location_flags = str2double(tab8(strcmp(tab8.Row,'flag'),:).Variables);
    success_criterion(end+1) = {location_flags};

    idx_remove = ismember(tab8.Row,{'flag','Optimal Material index'});
    tab8latex = tab8(~idx_remove,:);
    colspec = '>{\centering\arraybackslash}p{0.20\linewidth}>{\centering\arraybackslash}p{0.08\linewidth}>{\centering\arraybackslash}p{0.15\linewidth}>{\centering\arraybackslash}p{0.15\linewidth}>{\centering\arraybackslash}p{0.15\linewidth}>{\centering\arraybackslash}p{0.18\linewidth}';
    firstrow = '&& \multicolumn{4}{c}{Location}\\  \cline{3-6}';
    table2latex(tab8latex,[save_folder 'table_22.tex'],colspec,firstrow)
end


end
