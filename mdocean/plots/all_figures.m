% Generate all figures used in the report
function [num_figs, num_tabs] = all_figures( which_figs, which_tabs )

num_figs = 10;
num_tabs = 7;
if nargin==0
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

%% figure 1 - RM3 image
if any(which_figs == 1)
    % Created in powerpoint
end

%% figure 2 - N2 diagram
if any(which_figs == 2)
    % Created in powerpoint
end

%% figure 3, 4 - saturation time signal, saturation alpha
if any(which_figs == 3 | which_figs == 4)
    sin_saturation_demo()
end

%% figure 5 - JPD multiplication
if any(which_figs == 5)
    p = parameters();
    b = var_bounds(p);
    X = [b.X_noms; 1];
    plot_power_matrix(X,p)
end


%% figure 6, 7 - pareto front, design heuristics
if any(which_figs == 6 | which_figs == 7)  
    pareto_search();
    pareto_curve_heuristics()
end

%% figure 8 - parameter sensitivities
if any(which_figs == 8)
    param_sweep()
end

%% figure 9, 10 - overlaid geometry, probability CDF
%% table 5 - optimal DVs for 4 designs
if any(which_figs == 9 | which_figs == 10 | which_tabs == 5)
    compare()
end

%% table 1 - design variables table
if any(which_tabs == 1)
    p = parameters();
    b = var_bounds(p);
    DV_table = array2table([b.X_mins b.X_noms b.X_maxs], ...
        'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', b.var_names(1:end-1));
    display(DV_table)
end

%% table 2 - constraints table
if any(which_tabs == 2)
    % todo
end

%% table 3 - parameters table
if any(which_tabs == 3)
    p = parameters();
    struct2table(p, 'AsArray',true)
end

%% table 4 - validation table
if any(which_tabs == 4)
    [~,~,~,~,tab] = validate_nominal_RM3();
    display(tab)
end

%% table 6 - optimal DVs for 4 locations
if any(which_tabs == 6)
    location_sensitivity()
end

%% paragraph 4.2 - convergence for different x0s
if any(which_tabs == 7)
    %gradient_mult_x0()
end

end