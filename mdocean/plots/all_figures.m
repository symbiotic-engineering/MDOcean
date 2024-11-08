% Generate all figures used in the paper
function [success_criterion,num_figs,num_tabs,fig_names,tab_names] = all_figures( which_figs, which_tabs )

num_figs = 10;
num_tabs = 7;
fig_names = cell([1,num_figs]);
tab_names = cell([1,num_tabs]);
success_criterion = {};

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

%% figure 1 - RM3 image
fig_names{1} = 'Fig. 1: RM3 image';
if any(which_figs == 1)
    % Created in powerpoint
end

%% figure 2 - N2 diagram
fig_names{2} = 'Fig. 2: N2 diagram';
if any(which_figs == 2)
    % Created in powerpoint
end

%% figure 3, 4 - saturation time signal, saturation alpha
fig_names{3} = 'Fig. 3: saturation time signal';
fig_names{4} = 'Fig. 4: saturation alpha';
if any(which_figs == 3 | which_figs == 4)
    sin_saturation_demo()
end

%% figure 5 - JPD multiplication
fig_names{5} = 'Fig. 5: JPD multiplication';
if any(which_figs == 5)
    p = parameters();
    b = var_bounds();
    X = [b.X_noms; 1];
    plot_power_matrix(X,p)
end


%% figure 6, 7 - pareto front, design heuristics
fig_names{6} = 'Fig. 6: pareto front';
fig_names{7} = 'Fig. 7: design heuristics';
if any(which_figs == 6 | which_figs == 7)  
    pareto_search();
    pareto_curve_heuristics()
end

%% figure 8 - parameter sensitivities
fig_names{8} = 'Fig. 8: parameter sensitivities';
if any(which_figs == 8)
    param_sweep()
end

%% figure 9, 10 - overlaid geometry, probability CDF
fig_names{9} = 'Fig. 9: overlaid geometry';
fig_names{10} = 'Fig. 10: probability CDF';
if any(which_figs == 9 | which_figs == 10 | which_tabs == 5)
    tab_5 = compare();
end

%% table 1 - design variables table
tab_names{1} = 'Tab. 1: design variables';
if any(which_tabs == 1)
    b = var_bounds();
    DV_table = array2table([b.X_mins b.X_noms b.X_maxs], ...
        'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', b.var_names(1:end-1));
    display(DV_table)
end

%% table 2 - constraints table
tab_names{2} = 'Tab. 2: constraints';
if any(which_tabs == 2)
    display(b.constraint_names')
end

%% table 3 - parameters table
tab_names{3} = 'Tab. 3: parameters';
if any(which_tabs == 3)
    [~,tab] = parameters();
    display(tab)
end

%% table 4 - validation table
tab_names{4} = 'Tab. 4: validation against nominal';
if any(which_tabs == 4)
    [~,~,~,~,tab] = validate_nominal_RM3();
    display(tab)
end

%% table 5 - optimal DVs for 4 designs
tab_names{5} = 'Tab. 5: optimal DVs for 4 designs';
if any(which_tabs == 5)
    % computation above with figures 9-10
    display(tab_5);
end

%% table 6 - optimal DVs for 4 locations
tab_names{6} = 'Tab. 6: optimal DVs for 4 locations';
if any(which_tabs == 6)
    tab = location_sensitivity();
    display(tab);
    location_flags = tab(strcmp(tab.Row,'flag'),:).Variables;
    success_criterion(end+1) = {location_flags};
end

%% paragraph 4.2 - convergence for different x0s
tab_names{7} = 'Tab. 7: convergence for different x0s';
if any(which_tabs == 7)
    %gradient_mult_x0()
end

end