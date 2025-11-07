function fig = plot_power_matrix(X, p, b, filename_uuid)

    % get matrices to plot
    [CW_over_CW_max, P_wave, CW_max, P_elec, ...
        force_sat_ratio, drag_ratio, eff] = check_max_CW(filename_uuid, p, X, false);
    [T, Hs] = meshgrid(p.T, p.Hs);
    p.JPD(p.JPD > 0 & p.JPD < .001) = 0;
    hrs_in_yr = 8766;
    hours = p.JPD / 100 * hrs_in_yr;
    energy = P_elec .* hours;

    % check that multiplication is correct
    P_elec_calc = P_wave .* CW_max .* CW_over_CW_max .* drag_ratio .* force_sat_ratio .* eff;
    err = (P_elec - P_elec_calc) ./ P_elec;
    assert(all(abs(err(~isnan(err))) < 1e-3, 'all'));

    % figure
    fig = figure;
    t = tiledlayout(2, 12);
    t.TileSpacing = 'tight';
    t.Padding = 'compact';

    mycontour(T, Hs, P_wave / 1000, 'Wave Power (kW/m)', false, true);

    plot_char('x');

    mycontour(T, Hs, CW_max, 'Max Capture Width (m)', false, true);

    plot_char('x');

    mycontour(T, Hs, CW_over_CW_max * 100, 'Radiation Efficiency (%)', true, true);

    plot_char('x');

    mycontour(T, Hs, drag_ratio * 100, 'Drag Efficiency (%)', true, true);

    plot_char('x');

    mycontour(T, Hs, force_sat_ratio * 100, 'F_{max} Factor (%)', true, false);

    plot_char('x');

    mycontour(T, Hs, eff * 100, 'Electrical Efficiency (%)', true, false);

    plot_char('x');

    mycontour(T, Hs, p.JPD, 'Site Probability (%)', false, false);

    plot_char('=');

    mycontour(T, Hs, energy / 1e6, 'Annual Energy (MWh)', false, false);

    xlabel(t, 'Wave Period T_e (s)');
    ylabel(t, 'Wave Height H_s (m)');
    title(t, ' ');
    improvePlot;

    set(fig, 'Position', [0 123 1530 620]);

end

function mycontour(X, Y, Z, title_text, scale_100, title_higher)
    nexttile([1 2]);

    Z_is_constant = numel(unique(Z(~isnan(Z)))) == 1;
    if Z_is_constant
        % avoid "contour not rendered for constant zdata"
        x = [min(X, [], 'all') max(X, [], 'all')];
        y = [min(Y, [], 'all') max(Y, [], 'all')];
        imagesc('XData', x, 'YData', y, 'CData', Z, 'AlphaData', ~isnan(Z));
    else
        if scale_100
            levels = sort([min(Z(:)), max(Z(:)), 0:10:100]);
        else
            levels = 10;
        end
        contourf(X, Y, Z, levels);
    end
    if title_higher
        title(title_text, 'Position', [11.75, 7, 0]);
    else
        title(title_text);
    end
    colorbar;
    if scale_100
        caxis([0 100]);
    end
end

function plot_char(c)
    nexttile;
    text(0, 0.5, c, 'FontSize', 40);
    axis off;
end
