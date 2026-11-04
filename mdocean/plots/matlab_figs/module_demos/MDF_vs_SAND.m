function fig = MDF_vs_SAND()
    N_dyn_vec = [5 20 50 100 200];
    N_mdf_vec = [50 100 200 500 1000];

    t_dynam_over_t_other = logspace(-2,2);
    factor = linspace(0,1);
    [t_dynam_over_t_other_mesh, factor_mesh] = meshgrid(t_dynam_over_t_other, factor);

    fig = figure;
    t = tiledchartlayout(length(N_dyn_vec), length(N_mdf_vec));

    for i_Ndyn = 1:length(N_dyn_vec)
      for i_Nmdf = 1:length(N_mdf_vec)
        N_dyn = N_dyn_vec(i_Ndyn);
        N_mdf = N_mdf_vec(i_Nmdf);

        N_sand_min = min(N_dyn,N_mdf);
        N_sand_max = N_dyn + N_mdf;

        N_sand_mesh = N_sand_min + factor_mesh * (N_sand_max - N_sand_min);
        
        t_mdf_over_t_sand_mesh = N_mdf ./ N_sand_mesh * (N_dyn * t_dynam_over_t_other_mesh + 1) ./ (t_dynam_over_t_other_mesh + 1);

        t_dynam_over_t_total_mesh = 1./(1+ 1./t_dynam_over_t_other_mesh);
        nexttile;
        contourf(t_dynam_over_t_total_mesh, factor_mesh, t_mdf_over_t_sand_mesh)
        xlabel('t_{dynam}/t_{total}')
        ylabel('f_{N_{SAND}}')
        
    end

    xlabel(t, 'N_{dynam}')
    ylabel(t, 'N_{MDF}')

    title(t, 't_{MDF}/t_{SAND}')

    cb = colorbar;
    colormap(bluewhitered)
    cb.Layout.Tile = 'eastoutside';
end