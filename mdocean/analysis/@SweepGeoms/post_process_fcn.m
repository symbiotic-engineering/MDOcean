function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
% post_process_fcn  Generate all geometry-sweep figures.
%
% Produces three original figures (parallel-coordinates, scatter,
% line-plot), six capture-width-ratio (CWR) figures — scatter and
% line variants for each of three characteristic dimensions:
%   1. a_2      – float radius
%   2. V^{1/3}  – cube root of displaced volume  (val.vol_f + val.vol_s)
%   3. SA^{1/2} – square root of wetted surface area (cylindrical approx.)
% plus a Pareto plot (CW vs surface area), a grid scatter matrix
% (6 sweep variables × 8 output quantities), and a 6×6 factor-analysis
% pcolor of main and interaction η² effects on CW/SA^{1/2}.
%
% :param intermed_result_struct: Output of analysis_fcn.
% :returns: fig_array            12-element array of figure handles
% :returns: tab_array_display    Empty (no tables)
% :returns: tab_array_latex      Empty (no tables)
% :returns: end_result_struct    Struct with completion flag
% :returns: tab_firstrows        Empty
% :returns: tab_colspecs         Empty

    % ------------------------------------------------------------------
    % Unpack intermediate results
    % ------------------------------------------------------------------
    hydro_ratio_result = intermed_result_struct.hydro_ratio_result;
    m0h_stored         = intermed_result_struct.m0h_stored;
    val                = intermed_result_struct.val;
    H                  = intermed_result_struct.H;
    A1_A2              = intermed_result_struct.A1_A2;
    A2_H               = intermed_result_struct.A2_H;
    D1_H               = intermed_result_struct.D1_H;
    D2_D1              = intermed_result_struct.D2_D1;
    A3_A1              = intermed_result_struct.A3_A1;
    A1                 = intermed_result_struct.A1;
    A2                 = intermed_result_struct.A2;
    D1                 = intermed_result_struct.D1;
    D2                 = intermed_result_struct.D2;
    h_vec              = intermed_result_struct.h_vec;
    a2_h               = intermed_result_struct.a2_h;
    T                  = intermed_result_struct.p.T;
    g                  = intermed_result_struct.p.g;

    % ------------------------------------------------------------------
    % Compute CW (capture width) and CWR for three characteristic dims
    % ------------------------------------------------------------------
    nT = length(T);

    % CW_max(T) = g*T^2/(4*pi^2), broadcast over all geometry dimensions
    CW_max_vec  = g * T(:).^2 / (4*pi^2);                         % [nT x 1]
    CW_max_grid = reshape(CW_max_vec, [nT, ones(1, ndims(H))]);    % [nT x 1 x ...]
    CW          = hydro_ratio_result .* CW_max_grid;               % [nT x size(A1)]

    % 1. charac_dim = a_2 (float radius)
    charac_dim_a2  = reshape(A2, [1, size(A2)]);                   % [1 x size(A1)]
    CWR_a2         = CW ./ charac_dim_a2;

    % 2. charac_dim = cube root of displaced volume
    vol_f          = reshape([val.vol_f], size(A1));
    vol_s          = reshape([val.vol_s], size(A1));
    charac_dim_vol = reshape((vol_f + vol_s).^(1/3), [1, size(A1)]);
    CWR_vol        = CW ./ charac_dim_vol;

    % 3. charac_dim = SA^{1/2}: square root of wetted surface area (cylindrical approx.)
    %    Float: outer lateral cylinder + top annular cap
    %    Spar:  outer lateral cylinder + bottom disk cap
    D_f_grid = 2 * A2;
    D_s_grid = 2 * A1;
    SA_float      = pi * D_f_grid .* D2 + pi * (D_f_grid/2).^2;
    SA_spar       = pi * D_s_grid .* D1 + pi * (D_s_grid/2).^2;
    SA_total      = SA_float + SA_spar;
    charac_dim_sa = reshape(sqrt(SA_total), [1, size(A1)]);
    CWR_sa        = CW ./ charac_dim_sa;

    % ------------------------------------------------------------------
    % Shared variables for scatter / line plots
    % ------------------------------------------------------------------
    num_m0h = nT * length(h_vec);
    [m0h_mat, order]  = myreshape(m0h_stored,         num_m0h);
    result_mat        = myreshape(hydro_ratio_result,  num_m0h, order(:,1));
    CWR_a2_mat        = myreshape(CWR_a2,              num_m0h, order(:,1));
    CWR_vol_mat       = myreshape(CWR_vol,             num_m0h, order(:,1));
    CWR_sa_mat        = myreshape(CWR_sa,              num_m0h, order(:,1));

    % Colours / sizes (encode geometry dims as RGB)
    red       = myresize(A1_A2, nT);
    green     = myresize(D1_H,  nT);
    blue      = myresize(D2_D1, nT);
    color     = [red, green, blue];
    size_var  = myresize(A3_A1, nT);

    red2      = A1_A2(1,:,:,:,:,:);
    green2    = D1_H(1,:,:,:,:,:);
    blue2     = D2_D1(1,:,:,:,:,:);
    color2    = [red2(:), green2(:), blue2(:)];

    size_var2     = A3_A1(1,:,:,:,:,:);
    size_mult     = 3;
    size_var_name = 'a_3/a_1';

    marker_type_var = A2_H(1,:,:,:,:,:);
    marker_var_name = 'a_2/h';
    marker_types    = {'o','x','v','s','+','.'};

    m0h_minmax = [min(m0h_mat(:)), max(m0h_mat(:))];

    % ------------------------------------------------------------------
    % Figure 1: Parallel-coordinates plot
    % ------------------------------------------------------------------
    m0h_tmp    = m0h_stored(1,:,:,:,:,:,:);
    result_tmp = hydro_ratio_result(1,:,:,:,:,:,:);
    result_disc = discretize(result_tmp, 10, 'categorical');
    [result_disc_sorted, idx_sort] = sort(result_disc(:));
    idx = idx_sort(~isundefined(result_disc_sorted));

    % Use tiledlayout(1,1) so an overlay axes can be placed on the same tile
    % (same strategy as multistart_postpro.m lines 155-183) to apply LaTeX
    % tick labels — parallelplot has no LabelInterpreter property.
    fig1 = figure('Visible','off');
    t1   = tiledlayout(fig1, 1, 1, 'Padding', 'compact');
    nexttile(t1);
    pp = parallelplot([H(idx), A1_A2(idx), A2_H(idx), D1_H(idx), D2_D1(idx), ...
                       A3_A1(idx), m0h_tmp(idx), result_tmp(idx)], ...
                      'GroupData', result_disc_sorted(~isundefined(result_disc_sorted)));
    pp.CoordinateTickLabels = '';   % suppress original labels; overlay supplies them
    pp.Color = parula(10);

    coord_labels = {'$h$','$a_1/a_2$','$a_2/h$','$d_1/h$','$d_2/d_1$', ...
                    '$a_3/a_1$','$m_0 h$','$CW/CW_{\max}$'};
    nC = numel(coord_labels);
    hOvl = axes(t1);
    hOvl.Layout.Tile = 1;
    set(hOvl, 'XLim', [0.5, nC+0.5], 'Color', 'none', 'XColor', 'black', ...
              'YColor', 'none', 'XTick', 1:nC, 'XTickLabel', coord_labels, ...
              'TickLabelInterpreter', 'latex', 'Box', 'off', 'FontSize', 9);

    % ------------------------------------------------------------------
    % Figure 2: Scatter plot  (m0h vs CW/CW_max)
    % ------------------------------------------------------------------
    fig2 = figure('Visible','off');
    scatter(m0h_stored(:), hydro_ratio_result(:), size_var, color)
    set(gca, 'XScale', 'log')
    xlabel('m_0 h')
    ylabel('CW/CW_{max}')
    ylim([0 1])

    % ------------------------------------------------------------------
    % Figure 3: Line plot  (m0h vs CW/CW_max)
    % ------------------------------------------------------------------
    fig3 = make_line_fig(m0h_mat, result_mat, 'CW/CW_{max}', [0 1.5], ...
                         color2, size_var2, size_mult, size_var_name, ...
                         marker_type_var, marker_var_name, marker_types, a2_h, m0h_minmax);

    % ------------------------------------------------------------------
    % Figures 4-6: CWR scatter plots
    % ------------------------------------------------------------------
    fig4 = make_scatter_fig(m0h_stored, CWR_a2,  size_var, color, 'CW/a_2');
    fig5 = make_scatter_fig(m0h_stored, CWR_vol, size_var, color, 'CW/V^{1/3}');
    fig6 = make_scatter_fig(m0h_stored, CWR_sa,  size_var, color, 'CW/SA^{1/2}');

    % ------------------------------------------------------------------
    % Figures 7-9: CWR line plots
    % ------------------------------------------------------------------
    fig7 = make_line_fig(m0h_mat, CWR_a2_mat,  'CW/a_2',         [], ...
                         color2, size_var2, size_mult, size_var_name, ...
                         marker_type_var, marker_var_name, marker_types, a2_h, m0h_minmax);
    fig8 = make_line_fig(m0h_mat, CWR_vol_mat, 'CW/V^{1/3}',     [], ...
                         color2, size_var2, size_mult, size_var_name, ...
                         marker_type_var, marker_var_name, marker_types, a2_h, m0h_minmax);
    fig9 = make_line_fig(m0h_mat, CWR_sa_mat,  'CW/SA^{1/2}',   [], ...
                         color2, size_var2, size_mult, size_var_name, ...
                         marker_type_var, marker_var_name, marker_types, a2_h, m0h_minmax);

    % ------------------------------------------------------------------
    % Figure 10: Pareto plot  (CW vs Surface Area)
    % ------------------------------------------------------------------
    % Collapse the T dimension by taking the best CW for each geometry.
    CW_max_T = squeeze(max(CW, [], 1));   % [size(A1)]

    % Color/size convention extended to all geometry points (incl. both h
    % levels) so every point in the 486-geometry grid is represented.
    color_pareto    = [A1_A2(:), D1_H(:), D2_D1(:)];
    size_var_pareto = A3_A1(:);

    fig10 = make_pareto_fig(CW_max_T(:), SA_total(:), ...
                            color_pareto, size_var_pareto, size_mult);

    % ------------------------------------------------------------------
    % Figure 11: Grid scatter matrix
    %   Columns = 6 sweep vars (5 geometry ratios + m0h)
    %   Rows    = 8 output quantities
    % ------------------------------------------------------------------
    a1_a2_flat = myresize(A1_A2, nT);
    a2_h_flat  = myresize(A2_H,  nT);
    d1_h_flat  = myresize(D1_H,  nT);
    d2_d1_flat = myresize(D2_D1, nT);
    a3_a1_flat = myresize(A3_A1, nT);

    x_vars   = {a1_a2_flat, a2_h_flat, d1_h_flat, d2_d1_flat, a3_a1_flat, m0h_stored(:)};
    x_labels = {'a_1/a_2', 'a_2/h', 'd_1/h', 'd_2/d_1', 'a_3/a_1', 'm_0h'};

    charac_a2_flat  = myresize(A2,                       nT);
    charac_vol_flat = myresize((vol_f + vol_s).^(1/3),   nT);
    charac_sa_flat  = myresize(sqrt(SA_total),           nT);

    y_vars   = {hydro_ratio_result(:), charac_a2_flat, charac_vol_flat, charac_sa_flat, ...
                CW(:), CWR_a2(:), CWR_vol(:), CWR_sa(:)};
    y_labels = {'CW/CW_{max}', 'a_2 (m)', 'V^{1/3} (m)', 'SA^{1/2} (m)', ...
                'CW (m)', 'CW/a_2', 'CW/V^{1/3}', 'CW/SA^{1/2}'};

    fig11 = make_grid_scatter_fig(x_vars, x_labels, y_vars, y_labels, color, size_var);

    % ------------------------------------------------------------------
    % Figure 12: Factor analysis — 6×6 imagesc of η² effects on CW/SA^{1/2}
    %   Diagonal  (i==j): main effect of input i
    %   Off-diag  (i≠j): two-factor interaction effect between inputs i and j
    % ------------------------------------------------------------------
    Y_fa = CW_max_T ./ sqrt(SA_total);   % [size(A1)], per-geometry CW/SA^{1/2}
    factor_labels = {'h', 'a_1/a_2', 'a_2/h', 'd_1/h', 'd_2/d_1', 'a_3/a_1'};
    fig12 = make_factor_fig(Y_fa, factor_labels);

    fig_array = [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12];

    tab_array_display = {};
    tab_array_latex   = {};
    tab_firstrows     = {};
    tab_colspecs      = {};

    end_result_struct.sweep_geoms_complete = true;
end

% ======================================================================
% Local helper functions
% ======================================================================

function fig = make_pareto_fig(CW_vec, SA_vec, color_pareto, size_var_pareto, size_mult)
%MAKE_PARETO_FIG  Pareto front: capture width (x) vs surface area (y).
%   Finds the max-CW / min-SA Pareto front via an inline sweep (sort by CW
%   descending, keep any point that achieves a new minimum SA).
    CW_v = CW_vec(:);
    SA_v = SA_vec(:);
    [CW_sort, sort_idx] = sort(CW_v, 'descend');
    SA_sort = SA_v(sort_idx);
    min_SA  = Inf;
    pareto_mask = false(numel(CW_v), 1);
    for k = 1:numel(CW_sort)
        if CW_sort(k) > 0 && SA_sort(k) < min_SA
            pareto_mask(k) = true;
            min_SA = SA_sort(k);
        end
    end
    idxo = sort_idx(pareto_mask);

    fig = figure('Visible', 'off');
    scatter(CW_vec(:), SA_vec(:), size_mult * size_var_pareto(:), color_pareto, ...
            'HandleVisibility', 'off')
    hold on
    plot(CW_vec(idxo), SA_vec(idxo), 'ks', ...
         'MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1.5, ...
         'DisplayName', 'Pareto front')
    xlabel('CW (m)')
    ylabel('Surface Area (m^2)')
    legend('Location', 'best')
    improvePlot
    set(gca, 'XScale', 'log', 'YScale', 'log')
end

function fig = make_grid_scatter_fig(x_vars, x_labels, y_vars, y_labels, color, size_var)
%MAKE_GRID_SCATTER_FIG  Grid matrix of scatter plots.
%   Rows correspond to y-axis output quantities; columns to x-axis sweep
%   variables.  Uses the same colour/size encoding as the other sweep
%   figures.
    n_x = numel(x_vars);
    n_y = numel(y_vars);

    fig = figure('Visible', 'off');
    t = tiledlayout(n_y, n_x, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iy = 1:n_y
        for ix = 1:n_x
            ax = nexttile(t);
            scatter(ax, x_vars{ix}, y_vars{iy}, size_var, color)
            if iy == n_y
                xlabel(ax, x_labels{ix})
            else
                set(ax, 'XTickLabel', [])
            end
            if ix == 1
                ylabel(ax, y_labels{iy})
            else
                set(ax, 'YTickLabel', [])
            end
        end
    end

    fig.Position(3:4) = [1400, 900];
    set(fig, 'PaperPositionMode', 'auto');
end

function fig = make_scatter_fig(m0h_stored, CWR, size_var, color, ylabel_str)
%MAKE_SCATTER_FIG  Scatter plot of m0h vs a CWR quantity.
    fig = figure('Visible','off');
    scatter(m0h_stored(:), CWR(:), size_var, color)
    xlabel('m_0 h')
    ylabel(ylabel_str)
end

function fig = make_line_fig(m0h_mat, y_mat, ylabel_str, ylim_vals, ...
                              color2, size_var2, size_mult, size_var_name, ...
                              marker_type_var, marker_var_name, marker_types, ...
                              a2_h, m0h_minmax)
%MAKE_LINE_FIG  Semilog line plot of m0h vs y_mat with colour/marker legends.
    fig = figure('Visible','off');
    ax  = gca();

    % Set ColorOrder BEFORE plotting so each line gets its assigned RGB colour.
    ax.ColorOrder = color2;

    num_lines = size(y_mat, 2);
    h_data    = gobjects(num_lines, 1);

    for i = 1:num_lines
        idx_marker = marker_type_var(i) == a2_h;
        marker     = marker_types{idx_marker};
        h_data(i)  = semilogx(ax, m0h_mat(:,i), y_mat(:,i), ...
                               ['-' marker], 'MarkerSize', size_mult * size_var2(i));
        set(h_data(i).Annotation.LegendInformation, 'IconDisplayStyle', 'off')
        hold on
    end

    xlabel('m_0 h')
    ylabel(ylabel_str)
    if ~isempty(ylim_vals)
        ylim(ylim_vals)
    end
    % Reference line y=1; excluded from legend.
    plot(ax, m0h_minmax, [1 1], 'k--', 'HandleVisibility', 'off')

    % -- Legend 1: marker type encodes a_2/h (on main axes ax) --
    h_marker = gobjects(length(a2_h), 1);
    for i = 1:length(a2_h)
        h_marker(i) = plot(ax, NaN, NaN, ['k' marker_types{i}], ...
                           'DisplayName', num2str(a2_h(i)));
    end
    h_marker_leg = legend(ax, h_marker, 'AutoUpdate', 'off', 'Location', 'northwest');
    title(h_marker_leg, marker_var_name)

    % -- Legend 2: marker size encodes a_3/a_1 (on invisible overlay axes ah1) --
    ah1         = axes('position', get(ax, 'position'));
    unique_size = unique(size_var2);
    h_size      = gobjects(length(unique_size), 1);
    for i = 1:length(unique_size)
        % Plot on ah1 so this legend is associated with ah1, not ax.
        h_size(i) = plot(ah1, NaN, NaN, 'ko', ...
                         'MarkerSize', size_mult * unique_size(i), ...
                         'DisplayName', num2str(unique_size(i)));
        hold(ah1, 'on')
    end
    h_size_leg = legend(ah1, 'AutoUpdate', 'off', 'Location', 'northeast');
    title(h_size_leg, size_var_name)
    ah1.Visible = 'off';

    % -- Annotation: explain RGB colour encoding --
    annotation(fig, 'textbox', [0.28, 0.01, 0.55, 0.06], ...
               'String', 'Color: R = a_1/a_2,  G = d_1/h,  B = d_2/d_1', ...
               'EdgeColor', 'none', 'Interpreter', 'tex', 'FontSize', 9, ...
               'HorizontalAlignment', 'center')

    set([h_data; h_marker], 'MarkerFaceColor', 'none')
    fig.Position(3:4) = [1000 600];
    xlim(ax, m0h_minmax)
end

function [byfreq, order] = myreshape(A, len, order)
%MYRESHAPE  Reshape A to [len, numel(A)/len] and optionally sort rows.
    reshaped = reshape(A, len, []);
    if nargin < 3
        [byfreq, order] = sort(reshaped);
    else
        byfreq = reshaped(order, :);
    end
end

function collapsed = myresize(A, len)
%MYRESIZE  Repeat A 'len' times along a new trailing dimension, then
%   collapse everything into a column vector.
%   Used to add a period dimension onto geometric meshes.
    n_repeats = [ones(1, ndims(A)), len];
    repeated  = repmat(A, n_repeats);
    shifted   = shiftdim(repeated, ndims(A));
    collapsed = shifted(:);
end
