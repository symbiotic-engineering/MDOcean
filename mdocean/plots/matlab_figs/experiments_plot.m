function [figs, results_tab] = experiments_plot(b, X_ins, ratios, LCOE, cost, power, failed, pareto_results_struct)
    % EXPERIMENTS_PLOT Post processing for design space exploration

    if nargin < 8
        pareto_results_struct = [];
    end

    % create table for display
    results_tab = array2table(X_ins, 'VariableNames', b.var_names(1:end - 1));
    LCOE = LCOE';
    cost = cost';
    power = power';
    results_tab = addvars(results_tab, round(LCOE(LCOE ~= Inf), 2), round(cost(cost ~= Inf), 1), failed, ...
                          'NewVariableNames', {'LCOE ($/kWh)', 'c_v (%)', 'Failed Constraints'});
    disp(results_tab);

    % plot pareto curve for comparison, if pareto results exist
    if ~isempty(pareto_results_struct)
        pareto_figs = pareto_curve_heuristics(pareto_results_struct);
        pareto_fig_num = pareto_figs(3).Number;
        fig1 = figure(pareto_fig_num);
        plot(power / 1e3, cost, '*--');

        title('Design of Experiments Pareto Front');
        l = legend(b.var_names_pretty);
        improvePlot;
        l.Location = 'bestoutside';
    else
        fig1 = gobjects(1, 1);
    end
    %% sensitivities plot
    [ratios_sorted, idx] = sort(ratios);
    LCOE(1, :) = LCOE(1, 1); % fill in nominal LCOE results for each DV where it wasn't repeatedly tested
    cost(1, :) = cost(1, 1);

    fig2 = figure;
    t = tiledlayout(2, 1);
    t.TileSpacing = 'compact';

    % LCOE subplot
    ax1 = nexttile(1);
    cols = {'r:', 'r--', 'r-', 'r-.', 'r.', ...       % bulk dims
            'b:', 'b--', ...                       % PTO
            'g:', 'g--', 'g-', 'g-.', 'g.'};  % structural
    yline(LCOE(1, 1), 'LineWidth', 2, 'Color', 'k', 'HandleVisibility', 'off');
    hold on;
    for i = 1:size(cols, 2)
        temp_LCOE = LCOE(idx, :).';
        plot(ratios_sorted, temp_LCOE(i, :), cols{i});
        hold on;
    end

    ylab1 = ylabel('LCOE ($/kWh)');
    x_range = [1 / 3 3];
    axis(ax1, [x_range .6 1.1]);
    l = legend(b.var_names_pretty{1:end - 1});
    l.Location = 'northeastoutside';
    grid on;
    hold off;

    % cost subplot
    ax2 = nexttile(2);
    yline(cost(1, 1), 'LineWidth', 2, 'Color', 'k');
    hold on;
    for i = 1:size(cols, 2)
        temp_cost = cost(idx, :).';
        plot(ratios_sorted, temp_cost(i, :), cols{i});
    end

    ylab2 = ylabel('Structural & PTO Cost ($M)');
    grid on;

    % shared plot
    title(t, 'Design of Experiments Results', 'FontWeight', 'bold', 'FontSize', 20);
    grid on;
    linkaxes([ax1, ax2], 'x');
    xlabel('Design Variable Ratio (-)');
    xticklabels(ax1, {});
    xticks(ax2, xticks(ax1));
    improvePlot;
    ylab1.FontSize = 16.5;
    ylab2.FontSize = 16.5;
    xlim(x_range);
    fig2.Position(3:4) = [600  666]; % make taller

    figs = [fig1, fig2];
end
