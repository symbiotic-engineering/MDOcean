% Generate all figures used in the paper
function [success_criterion,fig_output,tab_output,...
            num_figs,num_tabs,fig_names,tab_names] = all_figures( which_figs, which_tabs, filename_uuid )

if nargin<3
    filename_uuid = ''; % required argument for anything running gradient_optim in parallel,
                        % since prob2struct needs unique filenames for code generation
end

num_figs = 29;
num_tabs = 7;
fig_names = cell([1,num_figs]);
tab_names = cell([1,num_tabs]);
success_criterion = {};

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

fig_output = gobjects(1, length(which_figs));

fig_output = cell(1, length(which_figs));

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

%% figure 1 - RM3 image
fig_names{1} = 'Fig. 1: RM3 image';
if any(which_figs == 1)
    % Created in powerpoint
    fig1 = figure;
    imshow(imread("geometry.png"),'Parent',axes(fig1));
    fig_output{which_figs==1} = fig1;
end

%% figure 2 - Methodology overview
fig_names{2} = 'Fig. 2: Methodology overview';
if any(which_figs == 2)
    % Created in powerpoint
    fig2 = figure;
    imshow(imread("methods_flowchart_2_cropped.jpg"),'Parent',axes(fig2));
    fig_output{which_figs==2} = fig2;
end

%% figure 3 - N2 diagram
fig_names{3} = 'Fig. 3: N2 diagram';
if any(which_figs == 3)
    % Created in powerpoint
    fig3 = figure;
    imshow(imread("simple_N2.jpg"),'Parent',axes(fig3));
    fig_output{which_figs==3} = fig3;
end

%% figure 4 - Dimensions
fig_names{4} = 'Fig. 4: Dimensions';
if any(which_figs == 4)
    % Created in powerpoint
    fig4 = figure;
    imshow(imread("dimensions.jpg"),'Parent',axes(fig4));
    fig_output{which_figs==4} = fig4;
end

%% figure 5 - MEEM geometry
fig_names{5} = 'Fig. 5: MEEM geometry';
if any(which_figs == 5)
    % Created in powerpoint
    fig5 = figure;
    imshow(imread("MEEM-dims-basic-2.jpg"),'Parent',axes(fig5));
    fig_output{which_figs==5} = fig5;
end

%% figure 6 - hydro coeffs vs freq
fig_names{6} = 'Fig. 6: hydro coeffs vs freq';
if any(which_figs == 6)
    fig6 = figure;
    % fixme - not implemented
    fig_output{which_figs==6} = fig6;
end

%% figure 7 - drag DF
fig_names{7} = 'Fig. 7: drag describing function';
if any(which_figs == 7)
    fig7 = figure;
    % fixme - not implemented
    fig_output{which_figs==7} = fig7;
end

%% figure 8, XX - saturation time signal, saturation alpha
fig_names{8} = 'Fig. 8: force saturation time signal';
%fig_names{4} = 'Fig. XX: saturation alpha';
if any(which_figs == 8) % | which_figs == 4)
    sin_desc_fcn_demo()
    figX = gcf;
    fig8 = figure(figX.Number-1);
    fig_output(which_figs==8) = fig8;
    %fig_output(which_figs==4) = figX;
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

%% figure 10 - WecSim error breakdown
fig_names{10} = 'Fig. 10: WECSim error breakdown';
if any(which_figs == 10)
    % fixme - not implemented
    fig10 = figure;
    fig_output(which_figs==10) = fig10;
end

%% figure 11 - FBD
fig_names{11} = 'Fig. 11: FBD';
if any(which_figs == 11)
    % Created in powerpoint
    fig11 = figure;
    imshow(imread("FBD.jpg"),'Parent',axes(fig11));
    fig_output{which_figs==11} = fig11;
end

%% figure 12 - cost vs N WEC
fig_names{12} = 'Fig. 12: cost vs N WEC';
if any(which_figs == 12)
    % fixme - not implemented
    fig12 = figure;
    fig_output(which_figs==12) = fig12;
end

%% figure 13 - sim runtime
fig_names{13} = 'Fig. 13: sim runtime';
if any(which_figs == 13)
    % fixme - not implemented
    fig13 = figure;
    fig_output(which_figs==13) = fig13;
end

%% figure 14 - Optimization process flowchart
fig_names{14} = 'Fig. 14: Optimization flowchart';
if any(which_figs == 14)
    % Created in powerpoint
    fig14 = figure;
    imshow(imread("optimization_process.jpg"),'Parent',axes(fig14));
    fig_output{which_figs==14} = fig14;
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
fig_names{22} = 'Fig. 22: pareto front';
fig_names{23} = 'Fig. 23: design heuristics';
fig_names{24} = 'Fig. 24: objective heuristics';
fig_names{25} = 'Fig. 25: pareto front with LCOE contours';
fig_names{26} = 'Fig. 26: constraint activity';
if any(which_figs == 22 | which_figs == 23 | which_figs == 24)  
    pareto_search(filename_uuid);
    pareto_curve_heuristics()
    fig24 = gcf;
    fig23 = figure(fig24.Number - 1);
    fig22 = figure(fig24.Number - 3);
    fig26 = figure(fig24.Number - 6); % constraint activity
    fig26.Position = [1 41 1536 844.8000];
    fig_output(which_figs==22) = fig22;
    fig_output(which_figs==23) = fig23;
    fig_output(which_figs==24) = fig24;
    % fixme: 25 not implemented
    fig_output(which_figs==26) = fig26;
end

%% figure 27, 28, 29 - overlaid geometry, hydro coeffs, probability CDF
fig_names{27} = 'Fig. 27: overlaid geometry';
fig_names{28} = 'Fig. 28: overlaid hydro coeffs';
fig_names{29} = 'Fig. 29: probability CDF';
if any(which_figs == 27 | which_figs == 28 | which_tabs == 29)
    tab5 = compare(filename_uuid);
    n = gcf().Number;
    fig29 = figure(n-1);
    fig27 = figure(n-2);
    fig_output(which_figs==9) = fig27;
    fig_output(which_figs==10) = fig29;
    % fixme: 28 not implemented
end

% fixme: 30 to 46 not implemented

%% table 1 - design variables table
tab_names{1} = 'Tab. 1: design variables';
if any(which_tabs == 1)
    b = var_bounds();
    tab1 = array2table([b.X_mins b.X_noms b.X_maxs], ...
        'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', b.var_names(1:end-1));
    display(tab1)
    tab_output{which_tabs==1} = tab1;
end

%% table 2 - constraints table
tab_names{2} = 'Tab. 2: constraints';
if any(which_tabs == 2)
    b = var_bounds();
    tab2 = b.constraint_names';
    display(tab2)
    tab_output{which_tabs==2} = tab2;
end

%% table 3 - parameters table
tab_names{3} = 'Tab. 3: parameters';
if any(which_tabs == 3)
    [~,tab3] = parameters();
    display(tab3)
    tab_output{which_tabs==3} = tab3;
end

%% table 4 - validation table
tab_names{4} = 'Tab. 4: validation against nominal';
if any(which_tabs == 4)
    [~,~,~,~,tab4a] = validate_nominal_RM3('report');
    display(tab4a)
    [~,~,~,~,tab4b] = validate_nominal_RM3('wecsim');
    display(tab4b)

    % merge table 4a and 4b while preserving row order
    sharedCols = intersect(tab4a.Properties.VariableNames, tab4b.Properties.VariableNames);
    tab4a.RowNum = (1:length(tab4a.Properties.RowNames))';
    tab4b.RowNum = length(tab4a.Properties.RowNames) + (1:length(tab4b.Properties.RowNames))';
    
    tab4 = outerjoin(tab4a, tab4b, 'Keys', [{'RowNum'},sharedCols], 'MergeKeys', true);
    tab4 = removevars(tab4,'RowNum');
    tab4.Properties.RowNames = [tab4a.Properties.RowNames; tab4b.Properties.RowNames];
    display(tab4)

    tab_output{which_tabs==4} = tab4;
end

%% table 5 - optimal DVs for 4 designs
tab_names{5} = 'Tab. 5: optimal DVs for 4 designs';
if any(which_tabs == 5)
    % computation above with figures 9-10
    display(tab5);
    tab_output{which_tabs==5} = tab5;
end

%% table 6 - optimal DVs for 4 locations
tab_names{6} = 'Tab. 6: optimal DVs for 4 locations';
if any(which_tabs == 6)
    tab6 = location_sensitivity(filename_uuid);
    display(tab6);
    tab_output{which_tabs==6} = tab6;
    location_flags = tab6(strcmp(tab6.Row,'flag'),:).Variables;
    success_criterion(end+1) = {location_flags};
end

%% paragraph 4.2 - convergence for different x0s
tab_names{7} = 'Tab. 7: convergence for different x0s';
if any(which_tabs == 7)
    tab7 = gradient_mult_x0(filename_uuid);
    tab_output{which_tabs==7} = tab7;
end

end