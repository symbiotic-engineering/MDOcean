function run_single(p,b,X)
% Run and plot a single design. Defaults to nominal design if no X input.

    if nargin==0
        clear;close all;clc
        p = parameters();
        b = var_bounds();
    end
    if nargin<3
        X = [b.X_noms; 1];
    end

    [LCOE, P_var, ~, g, val] = simulation(X,p)

    [feasible,~,failed] = is_feasible(g,X,p,b)

    plot_power_matrix(X,p,b,'')

    figure
    power_PDF(X,p)

    visualize_geometry(X,p)

end
