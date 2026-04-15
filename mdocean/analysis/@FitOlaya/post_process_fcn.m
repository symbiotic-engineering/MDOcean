function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
% Function post_process_fcn
%
% :param intermed_result_struct: Intermediate results struct with raw data and derived signals
% :returns: fig_array
% :returns: tab_array_display
% :returns: tab_array_latex
% :returns: end_result_struct
% :returns: tab_firstrows
% :returns: tab_colspecs

    raw             = intermed_result_struct.raw;
    case_2          = intermed_result_struct.case_2;
    case_4          = intermed_result_struct.case_4;
    idxs_alpha_cell = intermed_result_struct.idxs_alpha_cell;
    idxs_beta_cell  = intermed_result_struct.idxs_beta_cell;

    figs_original = plot_original_signals(raw);
    figs_case = plot_case_signals(case_2, case_4);

    slope_term = 'C*exp(logx).^(-m1)+A*exp(logx).^m3';
    x1_dip_term = 'abs(1-(exp(logx)/x1)).^m4';
    ft = fittype(['log( (' slope_term ') .* ' x1_dip_term ')'],...
        dependent='logy',independent='logx',coefficients=["C" "m1" "m3" "m4" "x1" "A"]);
    fo = fitoptions(ft);
    fo.Lower = [0 0 0 0 0 0];
    fo.StartPoint = [1.5, 1/5, 1/3, .5, 3, .5];

    figs_case4 = fit_and_plot_case4(case_4, case_2, idxs_alpha_cell, idxs_beta_cell, ft, fo);

    fig_array = [figs_original figs_case figs_case4];

    tab_array_display = {};
    tab_array_latex = {};

    tab_firstrows = {};
    tab_colspecs = {};

    end_result_struct.fit_olaya_analysis_complete = true;
end

function figs = plot_original_signals(raw)
    % Plot the original signals as published in the paper (figures 3b, 4b, 12a-12d)
    f1 = figure;
    plot(raw.kRb_3b_alpha_2, raw.f_3b_alpha_2, 'r')
    hold on
    plot(raw.kRb_3b_alpha_1, raw.f_3b_alpha_1, 'k')
    plot(raw.kRb_3b_alpha_pt51, raw.f_3b_alpha_pt51, 'b')
    xlabel('k R_b')
    ylabel('f bar')
    title('3b')

    f2 = figure;
    plot(raw.kRb_4b_alpha_8, raw.f_4b_alpha_8, 'r')
    hold on
    plot(raw.kRb_4b_alpha_2, raw.f_4b_alpha_2, 'b')
    plot(raw.kRb_4b_alpha_1, raw.f_4b_alpha_1, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('4b')

    f3 = figure;
    plot(raw.kRb_12a_alpha_1, raw.f_12a_alpha_1, 'r')
    hold on
    plot(raw.kRb_12a_alpha_pt75, raw.f_12a_alpha_pt75, 'b')
    plot(raw.kRb_12a_alpha_pt51, raw.f_12a_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12a')

    f4 = figure;
    plot(raw.kRb_12b_alpha_1, raw.f_12b_alpha_1, 'r')
    hold on
    plot(raw.kRb_12b_alpha_pt75, raw.f_12b_alpha_pt75, 'b')
    plot(raw.kRb_12b_alpha_pt51, raw.f_12b_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12b')

    f5 = figure;
    plot(raw.kRb_12c_alpha_1, raw.f_12c_alpha_1, 'r')
    hold on
    plot(raw.kRb_12c_alpha_pt75, raw.f_12c_alpha_pt75, 'b')
    plot(raw.kRb_12c_alpha_pt51, raw.f_12c_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12c')

    f6 = figure;
    plot(raw.kRb_12d_alpha_1, raw.f_12d_alpha_1, 'r')
    hold on
    plot(raw.kRb_12d_alpha_pt75, raw.f_12d_alpha_pt75, 'b')
    plot(raw.kRb_12d_alpha_pt51, raw.f_12d_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12d')

    figs = [f1 f2 f3 f4 f5 f6];
end

function figs = plot_case_signals(case_2, case_4)
    % Plot case 2 and case 4 signals with various normalizations

    % Check that 3b and 4b match for alpha = 1 and 2 after normalization
    f1 = figure;
    plot(case_2.kRb{3}, case_2.f{3}, 'r')
    hold on
    plot(case_2.kRb{2}, case_2.f{2}, 'k')
    plot(case_2.kRb{1}, case_2.f{1}, 'b')
    plot(case_2.kRb_4b{3}, case_2.f_4b_norm{3}, 'm*--')
    plot(case_2.kRb_4b{2}, case_2.f_4b_norm{2}, 'r*')
    plot(case_2.kRb_4b{1}, case_2.f_4b_norm{1}, 'k*')
    title('Comparison: 3b solid, 4b star *')
    xlabel('k R_b')
    ylabel('f / (\rho g \pi R_c^2)')
    ylim([0 10])

    % Plots for case 2
    f2 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha})
        hold on
    end
    ylim([0 10])
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 )')

    f3 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.alpha(i_alpha)^2)
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 \alpha^2)')

    f4 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}*case_2.alpha(i_alpha)^2, case_2.f{i_alpha} / case_2.alpha(i_alpha)^2)
        hold on
    end
    xlim([0 60])
    ylim([0 5])
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p \alpha^2')
    ylabel('f/(\rho g \pi R_c^2 \alpha^2)')

    f5 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.alpha(i_alpha)^2 ./ abs( besselh(0,case_2.kRx{i_alpha})))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 \alpha^2 H_0(k R_x))')

    f6 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x))')

    f7 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(-case_2.kh{i_alpha} * case_2.e1_over_h))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{-k e_1})')
    legend('location','best')

    f8 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(case_2.kh{i_alpha} * (case_2.e1_over_h-case_2.e2_over_h)))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{k (e_1-e_2)})')
    legend('location','best')

    f9 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(-case_2.kh{i_alpha} * case_2.e2_over_h))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{-k e_2})')
    legend('location','best')

    f10 = figure;
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}*case_2.alpha(i_alpha)^2, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(-case_2.kh{i_alpha} * case_2.e1_over_h))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlim([0 100])
    xlabel('kh R_x/R_p \alpha^2')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{-k e_1})')
    legend('location','best')

    % Plots for case 4
    f11 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            plot(case_4.kRb{i_alpha, j_beta}, case_4.f{i_alpha, j_beta})
            hold on
        end
        legend(string(num2cell(case_4.beta)))
        xlabel('k R_b')
        ylabel('f/(\rho g \pi R_c^2 )')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    f12 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            plot(case_4.kh_Rx_over_Rp{i_alpha, j_beta}, case_4.f{i_alpha, j_beta})
            hold on
        end
        legend(string(num2cell(case_4.beta)))
        xlabel('k h R_x/R_p')
        ylabel('f/(\rho g \pi R_c^2 )')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    f13 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            x = case_4.kh_Rx_over_Rp{i_alpha, j_beta};
            kh = case_4.kh{i_alpha, j_beta};
            N0 = .5 * (1 + sinh(2 * kh) ./ (2* kh));
            scaling = abs( besselh(0,case_4.kRx{i_alpha, j_beta})) .* cosh(kh) ./ sqrt(N0);
            plot(x, case_4.f{i_alpha, j_beta} .* scaling)
            hold on
        end
        legend(string(num2cell(case_4.beta)))
        xlabel('k h R_x/R_p')
        ylabel('f/(\rho g \pi R_c^2 ) * H_0(k R_x) cosh(kh)/sqrt(N_0) ')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    f14 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            x = case_4.kh_Rx_over_Rp{i_alpha, j_beta} / case_4.beta(j_beta);
            kh = case_4.kh{i_alpha, j_beta};
            N0 = .5 * (1 + sinh(2 * kh) ./ (2* kh));
            exponent = kh * case_4.e2_over_h;
            scaling = 1 ./ abs( besselh(0,case_4.kRx{i_alpha, j_beta})) ./ exp(-exponent) ./ sqrt(kh);
            y = case_4.f{i_alpha, j_beta} .* scaling;
            loglog(x, y)
            hold on
        end
        legend(string(num2cell(case_4.beta)),'location','best')
        xlabel('k h R_x/R_p / \beta')
        ylabel('$f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_2}\sqrt{kh})$ ','Interpreter','latex')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    f15 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            x = case_4.kh_Rx_over_Rp{i_alpha, j_beta} / case_4.beta(j_beta);
            kh = case_4.kh{i_alpha, j_beta};
            N0 = .5 * (1 + sinh(2 * kh) ./ (2* kh));
            exponent = kh * case_4.e1_over_h;
            scaling = 1 ./ abs( besselh(0,case_4.kRx{i_alpha, j_beta})) ./ exp(-exponent) ./ sqrt(kh);
            y = case_4.f{i_alpha, j_beta} .* scaling;
            loglog(x, y)
            hold on
        end
        legend(string(num2cell(case_4.beta)),'location','best')
        xlabel('k h R_x/R_p / \beta')
        ylabel('$f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}\sqrt{kh})$ ','Interpreter','latex')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    f16 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            x = case_4.kh_Rx_over_Rp{i_alpha, j_beta} / case_4.beta(j_beta);
            kh = case_4.kh{i_alpha, j_beta};
            N0 = .5 * (1 + sinh(2 * kh) ./ (2* kh));
            exponent = kh * case_4.e1_over_h;
            scaling = 1 ./ abs( besselh(0,case_4.kRx{i_alpha, j_beta})) ./ exp(-exponent);
            y = case_4.f{i_alpha, j_beta} .* scaling;
            loglog(x, y)
            hold on
        end
        legend(string(num2cell(case_4.beta)),'location','best')
        xlabel('k h R_x/R_p / \beta')
        ylabel('$f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1})$ ','Interpreter','latex')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    f17 = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            x = case_4.kh_Rx_over_Rp{i_alpha, j_beta} / case_4.beta(j_beta);
            kh = case_4.kh{i_alpha, j_beta};
            N0 = .5 * (1 + sinh(2 * kh) ./ (2* kh));
            exponent = kh * case_4.e1_over_h;
            scaling = case_4.beta(j_beta)^2 ./ abs( besselh(0,case_4.kRx{i_alpha, j_beta})) ./ exp(-exponent) ./ sqrt(kh);
            y = case_4.f{i_alpha, j_beta} .* scaling;
            loglog(x, y)
            hold on
        end
        legend(string(num2cell(case_4.beta)),'location','best')
        xlabel('k h R_x/R_p / \beta')
        ylabel('$f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}\sqrt{kh}/\beta^2)$ ','Interpreter','latex')
        title(['\alpha = ' num2str(case_4.alpha(i_alpha))])
    end

    figs = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17];
end

function figs = fit_and_plot_case4(case_4, case_2, idxs_alpha_cell, idxs_beta_cell, ft, fo)
    % Fit and plot case 4 data
    x_string = 'k h R_x/R_p \alpha^2 / \beta';
    y_string = '$f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}\sqrt{kh})$';

    figs_sep = plot_fit_and_coeffs_sep_alpha_beta(case_4.kh_Rx_over_Rp_alpha2_over_beta, case_4.f_over_H0_exp_sqrt_kh, ...
        case_4.alpha, case_4.beta, ft, fo, x_string, y_string);

    % Pre-compute base x and y for all-alpha combined plots
    fun = @(c, ia, ib) c / case_4.beta(ib);
    x_base = cellfun(fun, case_4.kh_Rx_over_Rp, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    fun = @(c, ia, ib) c ./ abs(besselh(0, case_4.kRx{ia,ib})) ./ exp(-case_4.kh{ia,ib} * case_4.e1_over_h) ./ case_4.kh{ia,ib} / (case_4.alpha(ia) * max(1, case_4.alpha(ia)));
    y_base = cellfun(fun, case_4.f, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    fun = @(c, ia, ib) compute_fit_base(c, case_4.alpha(ia), case_4.beta(ib), case_4.kh{ia,ib});
    y_pred_base = cellfun(fun, x_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    fun = @(c, ia, ib) c * case_4.alpha(ia)^4;
    x_f3 = cellfun(fun, x_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    % version 1: semilogx with alpha^2 x-scaling
    fun = @(c, ia, ib) c * case_4.alpha(ia)^2;
    x_f1_v1 = cellfun(fun, x_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    fun = @(c, ia, ib) c.^2 * case_4.beta(ib)^2;
    y_f1_v1 = cellfun(fun, y_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    % version 2: loglog with alpha^3 x-scaling
    fun = @(c, ia, ib) c * case_4.alpha(ia)^3;
    x_f1_v2 = cellfun(fun, x_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    fun = @(c, ia, ib) c.^2 * case_4.beta(ib) * (1 + case_4.beta(ib)/case_4.alpha(ia));
    y_f1_v2 = cellfun(fun, y_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    % Manual fit predictions in y_f1 spaces (for comparison overlay on combined plots)
    fun = @(c, ia, ib) c.^2 * case_4.beta(ib)^2;
    y_pred_f1_v1 = cellfun(fun, y_pred_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    fun = @(c, ia, ib) c.^2 * case_4.beta(ib) * (1 + case_4.beta(ib)/case_4.alpha(ia));
    y_pred_f1_v2 = cellfun(fun, y_pred_base, idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    x_label_v1 = 'k h R_x/R_p \alpha^2 / \beta';
    y_label_v1 = '$\beta^2 \left[f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}kh R_x/R_b \alpha)\right]^2$ ';
    x_label_v2 = 'k h R_x/R_p \alpha^3 / \beta';
    y_label_v2 = '$\beta(1+\beta/\alpha) \left[f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}kh R_x/R_b \alpha)\right]^2$ ';

    % version 1: all alpha combined (manual + auto fits with coeffs)
    [figs_v1, fits_v1] = plot_fit_and_coeffs_all_alpha_beta(x_f1_v1, y_f1_v1, x_f3, y_base, y_pred_base, ...
        case_4.alpha, case_4.beta, @loglog, x_label_v1, y_label_v1, [0 10], y_pred_f1_v1, ft, fo);

    % version 2: all alpha combined (manual + auto fits with coeffs)
    figs_v2 = plot_fit_and_coeffs_all_alpha_beta(x_f1_v2, y_f1_v2, x_f3, y_base, y_pred_base, ...
        case_4.alpha, case_4.beta, @loglog, x_label_v2, y_label_v2, [0 20], y_pred_f1_v2, ft, fo);

    % version 3: second-order fit in x^1.5 space
    x_new = x_f1_v2;
    y_new = y_f1_v2;

    wn = 3;
    z = 0.005;
    second_order_fit = @(xx) 2*sqrt((wn^1.5-xx.^1.5).^2 + (2*z*xx).^2) / wn^1.5;

    fun = @(y, x) y .* x.^1.5;
    y_new_x15 = cellfun(fun, y_new, x_new, 'UniformOutput', false);

    fun = @(x) second_order_fit(x);
    y_pred_new = cellfun(fun, x_new, 'UniformOutput', false);

    [~, figs_v3] = plot_fit_all_alpha_beta(x_new, y_new_x15, x_new, y_new_x15, y_pred_new, ...
        case_4.alpha, case_4.beta, @loglog, ...
        'k h R_x/R_p \alpha^3 / \beta', ...
        '$\beta(1+\beta/\alpha) \left[f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}kh R_x/R_b \alpha)\right]^2 x^{1.5}$ ', ...
        [], y_pred_new);

    % Overlay plot: case 4 v1 auto + case 2 data
    fig_overlay = plot_case2_case4_overlay_v1(case_4, case_2, x_f1_v1, y_f1_v1, ...
        fits_v1, x_label_v1, y_label_v1);

    figs = [figs_sep figs_v1 figs_v2 figs_v3 fig_overlay];
end

function figs = plot_fit_and_coeffs_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string)
    % Plots manual fit, auto fit, and fit coefficients vs alpha/beta for separate-alpha subplots.
    [~,    fig_manual] = plot_fit_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string, false);
    [fits, fig_auto]   = plot_fit_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string, true);
    figs_coeffs = plot_fit_coeffs_vs_alpha_beta(fits, alpha_vec, beta_vec);
    figs = [fig_manual fig_auto figs_coeffs];
end

function [figs, fits_auto] = plot_fit_and_coeffs_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
                                                    alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, ...
                                                    y_f1_label, y_f1_lim, y_f1_pred_cell, ft, fo)
    % Plots manual fit, auto fit, and fit coefficients vs alpha/beta for all-alpha combined plots.
    [~,    figs_manual] = plot_fit_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
                    alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell);
    [fits_auto, figs_auto] = plot_fit_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
                    alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell, ft, fo);
    figs_coeffs = plot_fit_coeffs_vs_alpha_beta(fits_auto, alpha_vec, beta_vec);
    figs = [figs_manual figs_auto figs_coeffs];
end

function figs = plot_fit_coeffs_vs_alpha_beta(fits,alpha_vec,beta_vec)
    cols = {'r','g','b','k'};
    % coeffs vs beta for each alpha
    coeff_names = coeffnames(fits{1,1});
    f1 = figure;
    for i_alpha=1:3
        for i_coeff = 1:length(coeff_names)
            subplot(1,length(coeff_names),i_coeff)
            coeff_name = coeff_names{i_coeff};
            coeff_vals = zeros(1,4);
            for j_beta = 1:4
                f = fits{i_alpha,j_beta};
                coeff_vals(j_beta) = f.(coeff_name);
            end
            alpha_str = ['\alpha=' num2str(alpha_vec(i_alpha))];
            plot(beta_vec,coeff_vals,cols{i_alpha},'DisplayName',alpha_str)
            xlabel('\beta')
            ylabel(coeff_name)
            hold on
        end
    end
    legend

    % coeffs vs alpha for each beta
    f2 = figure;
    for j_beta=1:4
        for i_coeff = 1:length(coeff_names)
            subplot(1,length(coeff_names),i_coeff)
            coeff_name = coeff_names{i_coeff};
            coeff_vals = zeros(1,3);
            for i_alpha = 1:3
                f = fits{i_alpha,j_beta};
                coeff_vals(i_alpha) = f.(coeff_name);
            end
            beta_str = ['\beta=' num2str(beta_vec(j_beta))];
            plot(alpha_vec,coeff_vals,cols{j_beta},'DisplayName',beta_str)
            hold on
            xlabel('\alpha')
            ylabel(coeff_name)
        end
    end
    legend

    figs = [f1 f2];
end

function [fits, fig] = plot_fit_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string, auto_fit_on)
    cols = {'r','g','b','k'};
    if auto_fit_on
        x_string = ['ln( ' x_string ')'];
        y_string = ['ln( ' y_string ')'];
    end
    fits = cell(size(alpha_vec));
    fig = figure;
    for i_alpha=1:3
        subplot(1,3,i_alpha)
        for j_beta = 1:4
            x = x_cell{i_alpha,j_beta};
            y = y_cell{i_alpha,j_beta};

            if ~auto_fit_on
                x_pred = logspace(-2,2);
                y_pred = perform_manual_fit(x_pred, alpha_vec(i_alpha), beta_vec(j_beta));
                h = loglog(x, y);
                hold on
                loglog(x_pred,abs(y_pred),'--','Color',h.Color,'HandleVisibility','off')
            else
                f = perform_auto_fit(x,y,ft,fo);
                h = plot(f,log(x),log(y));
                h(1).Color = cols{j_beta};
                h(2).Color = cols{j_beta};
                h(1).HandleVisibility = 'off';
                hold on
                fits{i_alpha,j_beta} = f;
            end
        end
        legend(string(num2cell(beta_vec)),'location','best')

        if i_alpha==2
            xlabel(x_string)
        else
            xlabel('')
        end
        if i_alpha==1
            ylabel(y_string,'Interpreter','latex')
        else
            ylabel('')
        end
        title(['\alpha = ' num2str(alpha_vec(i_alpha))] )
    end
end

function [fits, figs] = plot_fit_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
        alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell, ft, fo)
    % Optional args:
    %   y_f1_pred_cell (12th): manual-fit predictions in y_f1 space for overlay on f1
    %   ft, fo         (13th, 14th): fittype/fitoptions for auto fitting; returns fits when provided
    has_f1_pred  = nargin >= 12 && ~isempty(y_f1_pred_cell);
    do_auto_fit  = nargin >= 14 && ~isempty(fo);
    fits = cell(length(alpha_vec), length(beta_vec));

    % Distinct colors for beta and markers for alpha
    beta_colors  = {[0 0.447 0.741], [0.85 0.325 0.098], [0.466 0.674 0.188], ...
                    [0.494 0.184 0.557], [0.929 0.694 0.125], [0.301 0.745 0.933]};
    alpha_markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};

    f1 = figure; f2 = figure; f3 = figure;
    for i_alpha = 1:length(alpha_vec)
        for j_beta = 1:length(beta_vec)
            alpha_beta_label = ['\alpha=' num2str(alpha_vec(i_alpha)) ', \beta=' num2str(beta_vec(j_beta))];
            x_f1 = x_f1_cell{i_alpha, j_beta};
            y_f1 = y_f1_cell{i_alpha, j_beta};

            col = beta_colors{mod(j_beta-1, length(beta_colors)) + 1};
            mkr = alpha_markers{mod(i_alpha-1, length(alpha_markers)) + 1};

            figure(f1)
            if has_f1_pred
                % data as markers WITH legend entry; manual fit as dashed (no legend entry)
                f1_plot_fcn(x_f1, y_f1, mkr, 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'DisplayName', alpha_beta_label);
                hold on
                f1_plot_fcn(x_f1, y_f1_pred_cell{i_alpha, j_beta}, '--', ...
                    'Color', col, 'HandleVisibility', 'off')
            else
                f1_plot_fcn(x_f1, y_f1, ['-' mkr], 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'DisplayName', alpha_beta_label)
            end
            hold on

            if do_auto_fit
                f = perform_auto_fit(x_f1, y_f1, ft, fo);
                fits{i_alpha, j_beta} = f;
                % overlay auto fit as solid line; legend entry is the data marker above
                x_fit_line = exp(linspace(log(min(x_f1(x_f1>0))), log(max(x_f1)), 100));
                y_fit_line = exp(feval(f, log(x_fit_line)));
                f1_plot_fcn(x_fit_line, y_fit_line, '-', 'Color', col, ...
                    'HandleVisibility', 'off');
                y_pred_ij = exp(feval(f, log(x_f1)));
            else
                y_pred_ij = y_pred_cell{i_alpha, j_beta};
            end

            figure(f2)
            scatter(y_cell{i_alpha, j_beta}, y_pred_ij, 36, col, mkr, 'DisplayName', alpha_beta_label)
            hold on

            figure(f3)
            scatter(x_f3_cell{i_alpha, j_beta}, y_pred_ij ./ y_cell{i_alpha, j_beta} - 1, 36, col, mkr, 'DisplayName', alpha_beta_label)
            hold on
        end
    end

    figure(f1)
    legend('location','eastoutside')
    xlabel(x_f1_label, 'FontSize', 14)
    ylabel(y_f1_label, 'Interpreter', 'latex', 'FontSize', 14)
    title('All case 4 data with fits')
    if ~isempty(y_f1_lim)
        ylim(y_f1_lim)
    end
    improvePlot

    figure(f2)
    plot([.001 10],[.001 10],'k','LineWidth',2,'DisplayName','Equality')
    legend('location','best')
    xlabel('Actual')
    ylabel('Predicted ','Interpreter','latex')
    title('Error')
    ylim([0 13])
    ax = gca();
    ax.XScale = 'log';
    ax.YScale = 'log';
    legend("Position",[0.74964,0.16095,0.23964,0.47976])

    figure(f3)
    legend('location','best')
    xlabel('k h R_x/R_p \alpha^2 / \beta')
    ylabel('predicted/actual - 1','Interpreter','latex')
    title('Fractional Error')
    ylim([0 2.5])

    figs = [f1 f2 f3];
end

function fig = plot_case2_case4_overlay_v1(case_4, case_2, x_f1_v1, y_f1_v1, ...
        fits_v1, x_label_v1, y_label_v1)
    % Overlay plot: duplicate the case 4 v1 auto plot and add case 2 data

    alpha_markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};
    beta_colors  = {[0 0.447 0.741], [0.85 0.325 0.098], [0.466 0.674 0.188], ...
                    [0.494 0.184 0.557], [0.929 0.694 0.125], [0.301 0.745 0.933]};

    fig = figure;

    % Re-plot case 4 data and auto fits (duplicating the v1 auto plot)
    for i_alpha = 1:length(case_4.alpha)
        mkr = alpha_markers{mod(i_alpha-1, length(alpha_markers)) + 1};
        for j_beta = 1:length(case_4.beta)
            col = beta_colors{mod(j_beta-1, length(beta_colors)) + 1};
            x_f1 = x_f1_v1{i_alpha, j_beta};
            y_f1 = y_f1_v1{i_alpha, j_beta};

            % Case 4 data as filled markers (legend entry)
            case4_label = sprintf('Case 4 (h/R_b=%g): \\alpha=%g, \\beta=%g', ...
                case_4.h_over_Rb, case_4.alpha(i_alpha), case_4.beta(j_beta));
            semilogx(x_f1, y_f1, mkr, 'Color', col, 'MarkerFaceColor', col, ...
                'MarkerSize', 4, 'DisplayName', case4_label);
            hold on

            % Case 4 auto fit lines (no legend entry)
            f = fits_v1{i_alpha, j_beta};
            x_fit_line = exp(linspace(log(min(x_f1(x_f1>0))), log(max(x_f1)), 100));
            y_fit_line = exp(feval(f, log(x_fit_line)));
            semilogx(x_fit_line, y_fit_line, '-', 'Color', col, 'HandleVisibility', 'off');
        end
    end

    % Overlay case 2 data transformed to v1 space
    % Case 2 has beta=1 effectively (normalized from fig 3b/4b)
    beta_c2 = 1;
    case2_marker_size = 8;
    for i_alpha = 1:length(case_2.alpha)
        alpha_c2 = case_2.alpha(i_alpha);
        kh_c2 = case_2.kh{i_alpha};
        kRx_c2 = case_2.kRx{i_alpha};

        % x_v1 = kh_Rx_over_Rp * alpha^2 / beta
        x_base_c2 = case_2.kh_Rx_over_Rp{i_alpha} / beta_c2;
        x_v1_c2 = x_base_c2 * alpha_c2^2;

        % y_base = f / (|H0(kRx)| * exp(-kh*e1/h) * kh * alpha * max(1,alpha))
        y_base_c2 = case_2.f{i_alpha} ./ abs(besselh(0, kRx_c2)) ...
            ./ exp(-kh_c2 * case_2.e1_over_h) ./ kh_c2 ...
            / (alpha_c2 * max(1, alpha_c2));
        % y_v1 = y_base^2 * beta^2
        y_v1_c2 = y_base_c2.^2 * beta_c2^2;

        % Use same alpha marker; color comes from beta (beta_c2=1, use distinct gray)
        mkr_c2 = alpha_markers{mod(i_alpha-1, length(alpha_markers)) + 1};
        col_c2 = [0.5 0.5 0.5]; % gray for case 2 (beta=1, not in case 4 beta set)
        case2_label = sprintf('Case 2 (h/R_b=%g): \\alpha=%g, \\beta=%g', ...
            case_2.h_over_Rb, alpha_c2, beta_c2);
        semilogx(x_v1_c2, y_v1_c2, mkr_c2, 'Color', col_c2, 'MarkerFaceColor', col_c2, ...
            'MarkerSize', case2_marker_size, ...
            'DisplayName', case2_label);
    end

    legend('location', 'eastoutside')
    xlabel(x_label_v1, 'FontSize', 14)
    ylabel(y_label_v1, 'Interpreter', 'latex', 'FontSize', 14)
    title('Case 4 fits with Case 2 data overlay')
    ylim([0 10])
    improvePlot
end

function f = perform_auto_fit(x,y,ft,fo)
    y(y<=0) = 1e-5;
    not_nan = ~isnan(x) & ~isnan(y) & x>0;
    f = fit(log(x(not_nan)), log(y(not_nan)), ft, fo);
end

function [x1, A0] = get_fit_params(alpha, beta)
    if alpha==1
        x1 = 2.7;
        A0 = .2;
    elseif alpha==0.75
        x1 = 10;
        if beta==0.5
            A0 = 0;
        else
            A0 = 0.03;
        end
    elseif alpha==0.51
        A0 = 0;
        x1 = 80;
    end
end

function y_pred = perform_manual_fit(x_pred,alpha,beta)
    [x1_0, A0] = get_fit_params(alpha, beta);
    x1 = x1_0*alpha^2;
    y_pred = piecewise_power_law(x_pred, x1, A0, beta);
end

function y_pred = compute_fit_base(x, alpha, beta, kh)
    [x1, A0] = get_fit_params(alpha, beta);
    y_pred = piecewise_power_law(x, x1, A0, beta) ./ sqrt(kh) / (alpha * max(1, alpha));
end

function y = piecewise_power_law(x, x1, A0, beta)
    C = 1.5; % left intercept
    m1 = .2; % left slope
    m2 = 1.8;%4.5*case_4_beta(j_beta);
    m3 = .3; % right slope

    A = A0/beta^2;

    left_slope_term = C*x.^-m1;
    right_slope_term = A*x.^m3;

    left_dip_term = (1-(x/x1).^m2).^2;
    right_dip_term = (1-(x1./x).^m2).^2;

    switch_term = double(x > x1);

    y = left_slope_term .* (1-switch_term) .* left_dip_term + ...
        right_slope_term .* switch_term .* right_dip_term;
end
