function figs = run_single(p,b,X)
% Run and plot a single design. Defaults to nominal design if no X input.

% Returns:
%  figs - array of figure handles created by this routine

    if nargin==0
        clear;close all;clc
        p = parameters();
        b = var_bounds();
    end
    if nargin<3
        X = [b.X_noms; 1];
    end

    [J, ~, g, val] = simulation(X,p)

    LCOE = J(1);
    P_var = J(2);

    [feasible,~,failed] = is_feasible(g,X,p,b)

    h_power = plot_power_matrix(X,p,b,'');

    h_pdf = figure;
    power_PDF(X,p)

    h_geom = visualize_geometry(X,p);

    figs = [h_power, h_pdf, h_geom];

end
