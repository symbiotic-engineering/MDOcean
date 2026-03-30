clear all
% First just replicate graphs in paper

% Load raw digitized data from CSV files
raw = load_raw_data();

% Plot original signals as published in the paper
plot_original_signals(raw);

% Create derived signals from raw data
[case_2, case_4, idxs_alpha_cell, idxs_beta_cell] = create_signals(raw);

% Plot case signals with various normalizations
plot_case_signals(case_2, case_4);

% Define fit type and options for case 4
slope_term = 'C*exp(logx).^(-m1)+A*exp(logx).^m3';
x1_dip_term = 'abs(1-(exp(logx)/x1)).^m4';
ft = fittype(['log( (' slope_term ') .* ' x1_dip_term ')'],...
    dependent='logy',independent='logx',coefficients=["C" "m1" "m3" "m4" "x1" "A"]);
fo = fitoptions(ft);
fo.Lower = [0 0 0 0 0 0];
fo.StartPoint = [1.5, 1/5, 1/3, .5, 3, .5];

% Fit and plot case 4
fit_and_plot_case4(case_4, idxs_alpha_cell, idxs_beta_cell, ft, fo);


% functions
function raw = load_raw_data()
    % Load raw digitized data from CSV files for figures 3b, 4b, 12a-12d
    data_3b_alpha_2    = readmatrix("3b_red");
    data_3b_alpha_1    = readmatrix("3b_black");
    data_3b_alpha_pt51 = readmatrix("3b_blue");

    raw.kRb_3b_alpha_2    = data_3b_alpha_2(:,1);    raw.f_3b_alpha_2    = data_3b_alpha_2(:,2);
    raw.kRb_3b_alpha_1    = data_3b_alpha_1(:,1);    raw.f_3b_alpha_1    = data_3b_alpha_1(:,2);
    raw.kRb_3b_alpha_pt51 = data_3b_alpha_pt51(:,1); raw.f_3b_alpha_pt51 = data_3b_alpha_pt51(:,2);

    data_4b_alpha_8 = readmatrix("4b_red");
    data_4b_alpha_2 = readmatrix("4b_blue");
    data_4b_alpha_1 = readmatrix("4b_black");

    raw.kRb_4b_alpha_8 = data_4b_alpha_8(:,1); raw.f_4b_alpha_8 = data_4b_alpha_8(:,2);
    raw.kRb_4b_alpha_2 = data_4b_alpha_2(:,1); raw.f_4b_alpha_2 = data_4b_alpha_2(:,2);
    raw.kRb_4b_alpha_1 = data_4b_alpha_1(:,1); raw.f_4b_alpha_1 = data_4b_alpha_1(:,2);

    data_12a_alpha_1    = readmatrix("12a_red");
    data_12a_alpha_pt75 = readmatrix("12a_blue");
    data_12a_alpha_pt51 = readmatrix("12a_black");

    raw.kRb_12a_alpha_1    = data_12a_alpha_1(:,1);    raw.f_12a_alpha_1    = data_12a_alpha_1(:,2);
    raw.kRb_12a_alpha_pt75 = data_12a_alpha_pt75(:,1); raw.f_12a_alpha_pt75 = data_12a_alpha_pt75(:,2);
    raw.kRb_12a_alpha_pt51 = data_12a_alpha_pt51(:,1); raw.f_12a_alpha_pt51 = data_12a_alpha_pt51(:,2);

    data_12b_alpha_1    = readmatrix("12b_red");
    data_12b_alpha_pt75 = readmatrix("12b_blue");
    data_12b_alpha_pt51 = readmatrix("12b_black");

    raw.kRb_12b_alpha_1    = data_12b_alpha_1(:,1);    raw.f_12b_alpha_1    = data_12b_alpha_1(:,2);
    raw.kRb_12b_alpha_pt75 = data_12b_alpha_pt75(:,1); raw.f_12b_alpha_pt75 = data_12b_alpha_pt75(:,2);
    raw.kRb_12b_alpha_pt51 = data_12b_alpha_pt51(:,1); raw.f_12b_alpha_pt51 = data_12b_alpha_pt51(:,2);
    raw.f_12b_alpha_pt75(raw.f_12b_alpha_pt75<0) = 1e-3;

    data_12c_alpha_1    = readmatrix("12c_red");
    data_12c_alpha_pt75 = readmatrix("12c_blue");
    data_12c_alpha_pt51 = readmatrix("12c_black");

    raw.kRb_12c_alpha_1    = data_12c_alpha_1(:,1);    raw.f_12c_alpha_1    = data_12c_alpha_1(:,2);
    raw.kRb_12c_alpha_pt75 = data_12c_alpha_pt75(:,1); raw.f_12c_alpha_pt75 = data_12c_alpha_pt75(:,2);
    raw.kRb_12c_alpha_pt51 = data_12c_alpha_pt51(:,1); raw.f_12c_alpha_pt51 = data_12c_alpha_pt51(:,2);

    data_12d_alpha_1    = readmatrix("12d_red");
    data_12d_alpha_pt75 = readmatrix("12d_blue");
    data_12d_alpha_pt51 = readmatrix("12d_black");

    raw.kRb_12d_alpha_1    = data_12d_alpha_1(:,1);    raw.f_12d_alpha_1    = data_12d_alpha_1(:,2);
    raw.kRb_12d_alpha_pt75 = data_12d_alpha_pt75(:,1); raw.f_12d_alpha_pt75 = data_12d_alpha_pt75(:,2);
    raw.kRb_12d_alpha_pt51 = data_12d_alpha_pt51(:,1); raw.f_12d_alpha_pt51 = data_12d_alpha_pt51(:,2);
    raw.f_12d_alpha_pt51(raw.f_12d_alpha_pt51<0) = 1e-3;
end

function plot_original_signals(raw)
    % Plot the original signals as published in the paper (figures 3b, 4b, 12a-12d)
    figure
    plot(raw.kRb_3b_alpha_2, raw.f_3b_alpha_2, 'r')
    hold on
    plot(raw.kRb_3b_alpha_1, raw.f_3b_alpha_1, 'k')
    plot(raw.kRb_3b_alpha_pt51, raw.f_3b_alpha_pt51, 'b')
    xlabel('k R_b')
    ylabel('f bar')
    title('3b')

    figure
    plot(raw.kRb_4b_alpha_8, raw.f_4b_alpha_8, 'r')
    hold on
    plot(raw.kRb_4b_alpha_2, raw.f_4b_alpha_2, 'b')
    plot(raw.kRb_4b_alpha_1, raw.f_4b_alpha_1, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('4b')

    figure
    plot(raw.kRb_12a_alpha_1, raw.f_12a_alpha_1, 'r')
    hold on
    plot(raw.kRb_12a_alpha_pt75, raw.f_12a_alpha_pt75, 'b')
    plot(raw.kRb_12a_alpha_pt51, raw.f_12a_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12a')

    figure
    plot(raw.kRb_12b_alpha_1, raw.f_12b_alpha_1, 'r')
    hold on
    plot(raw.kRb_12b_alpha_pt75, raw.f_12b_alpha_pt75, 'b')
    plot(raw.kRb_12b_alpha_pt51, raw.f_12b_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12b')

    figure
    plot(raw.kRb_12c_alpha_1, raw.f_12c_alpha_1, 'r')
    hold on
    plot(raw.kRb_12c_alpha_pt75, raw.f_12c_alpha_pt75, 'b')
    plot(raw.kRb_12c_alpha_pt51, raw.f_12c_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12c')

    figure
    plot(raw.kRb_12d_alpha_1, raw.f_12d_alpha_1, 'r')
    hold on
    plot(raw.kRb_12d_alpha_pt75, raw.f_12d_alpha_pt75, 'b')
    plot(raw.kRb_12d_alpha_pt51, raw.f_12d_alpha_pt51, 'k')
    xlabel('k R_b')
    ylabel('f bar')
    title('12d')
end

function [case_2, case_4, idxs_alpha_cell, idxs_beta_cell] = create_signals(raw)
    % Create case datasets and derived signals from raw digitized data
    % fig 3 and 4: case 2; fig 12: case 4
    case_2.h_over_Rb = 1/.2;
    case_2.e1_over_h = 0.25;
    case_2.e2_over_h = 0.35;
    case_4.h_over_Rb = 1;
    case_4.e1_over_h = 0.4;
    case_4.e2_over_h = 0.5;

    % Normalize 4b data by Rc^2 to match 3b convention
    beta_4b = 0.5;
    case_2.kRb_4b    = {raw.kRb_4b_alpha_1, raw.kRb_4b_alpha_2, raw.kRb_4b_alpha_8};
    case_2.f_4b_norm = {raw.f_4b_alpha_1 * 1^2 / beta_4b^2, ...
                        raw.f_4b_alpha_2 * 2^2 / beta_4b^2, ...
                        raw.f_4b_alpha_8 * 8^2 / beta_4b^2};

    % Consolidate 3b and 4b into one dataset (case 2) and
    % 12a, 12b, 12c, 12d into another dataset (case 4)
    case_2.alpha = [.51 1 2 8];
    case_2.kRb = {raw.kRb_3b_alpha_pt51  raw.kRb_3b_alpha_1  raw.kRb_3b_alpha_2  raw.kRb_4b_alpha_8};
    case_2.f   = {raw.f_3b_alpha_pt51    raw.f_3b_alpha_1    raw.f_3b_alpha_2    case_2.f_4b_norm{3}};

    case_4.alpha = [.51 .75 1];
    case_4.beta  = [.5 .35 .25 .15];
    % rows: alpha, cols: beta
    case_4.kRb = {raw.kRb_12a_alpha_pt51  raw.kRb_12b_alpha_pt51  raw.kRb_12c_alpha_pt51  raw.kRb_12d_alpha_pt51; ...
                  raw.kRb_12a_alpha_pt75  raw.kRb_12b_alpha_pt75  raw.kRb_12c_alpha_pt75  raw.kRb_12d_alpha_pt75; ...
                  raw.kRb_12a_alpha_1     raw.kRb_12b_alpha_1     raw.kRb_12c_alpha_1     raw.kRb_12d_alpha_1};
    case_4.f   = {raw.f_12a_alpha_pt51    raw.f_12b_alpha_pt51    raw.f_12c_alpha_pt51    raw.f_12d_alpha_pt51; ...
                  raw.f_12a_alpha_pt75    raw.f_12b_alpha_pt75    raw.f_12c_alpha_pt75    raw.f_12d_alpha_pt75; ...
                  raw.f_12a_alpha_1       raw.f_12b_alpha_1       raw.f_12c_alpha_1       raw.f_12d_alpha_1};

    % Index cells for cellfun (row = alpha index, col = beta index)
    idxs_alpha_cell = num2cell([1 1 1 1; 2 2 2 2; 3 3 3 3]);
    idxs_beta_cell  = num2cell([1 2 3 4; 1 2 3 4; 1 2 3 4]);

    % Create k h R_x / R_p = kRb * h/Rb * Rx/Rp (used as x axis)
    case_2.Rx_over_Rp = max(1./case_2.alpha, 1);
    fun = @(c, ia) c .* case_2.Rx_over_Rp(ia) * case_2.h_over_Rb;
    case_2.kh_Rx_over_Rp = cellfun(fun, case_2.kRb, num2cell(1:4), 'UniformOutput', false);

    case_4.Rx_over_Rp = max(1./case_4.alpha, 1);
    fun = @(c, ia) c .* case_4.Rx_over_Rp(ia) * case_4.h_over_Rb;
    case_4.kh_Rx_over_Rp = cellfun(fun, case_4.kRb, idxs_alpha_cell, 'UniformOutput', false);

    % Create k R_x
    case_2.Rx_over_Rb = max(case_2.alpha, 1);
    fun = @(c, ia) c .* case_2.Rx_over_Rb(ia);
    case_2.kRx = cellfun(fun, case_2.kRb, num2cell(1:4), 'UniformOutput', false);

    case_4.Rx_over_Rb = max(case_4.alpha, 1);
    fun = @(c, ia) c .* case_4.Rx_over_Rb(ia);
    case_4.kRx = cellfun(fun, case_4.kRb, idxs_alpha_cell, 'UniformOutput', false);

    % Create kh
    fun = @(c) c .* case_2.h_over_Rb;
    case_2.kh = cellfun(fun, case_2.kRb, 'UniformOutput', false);

    fun = @(c) c .* case_4.h_over_Rb;
    case_4.kh = cellfun(fun, case_4.kRb, 'UniformOutput', false);

    % Create k h R_x / R_p * alpha^2 / beta
    fun = @(c, ia, ib) c .* case_4.alpha(ia) ./ case_4.beta(ib);
    case_4.kh_Rx_over_Rp_alpha2_over_beta = cellfun(fun, case_4.kh_Rx_over_Rp, ...
        idxs_alpha_cell, idxs_beta_cell, 'UniformOutput', false);

    % Create f / ( |H_0(k R_x)| exp(-k e_1) sqrt(kh) )
    fun = @(c, ia, ib) c ./ ( abs( besselh(0, case_4.kRx{ia, ib}) ) .* ...
                               exp(-case_4.kh{ia, ib} * case_4.e1_over_h) .* ...
                               sqrt(case_4.kh{ia, ib}) );
    case_4.f_over_H0_exp_sqrt_kh = cellfun(fun, case_4.f, idxs_alpha_cell, idxs_beta_cell, ...
        'UniformOutput', false);
end

function plot_case_signals(case_2, case_4)
    % Plot case 2 and case 4 signals with various normalizations

    % Check that 3b and 4b match for alpha = 1 and 2 after normalization
    figure
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
    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha})
        hold on
    end
    ylim([0 10])
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 )')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.alpha(i_alpha)^2)
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 \alpha^2)')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}*case_2.alpha(i_alpha)^2, case_2.f{i_alpha} / case_2.alpha(i_alpha)^2)
        hold on
    end
    xlim([0 60])
    ylim([0 5])
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p \alpha^2')
    ylabel('f/(\rho g \pi R_c^2 \alpha^2)')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.alpha(i_alpha)^2 ./ abs( besselh(0,case_2.kRx{i_alpha})))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 \alpha^2 H_0(k R_x))')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x))')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(-case_2.kh{i_alpha} * case_2.e1_over_h))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{-k e_1})')
    legend('location','best')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(case_2.kh{i_alpha} * (case_2.e1_over_h-case_2.e2_over_h)))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{k (e_1-e_2)})')
    legend('location','best')

    figure
    for i_alpha=1:4
        plot(case_2.kh_Rx_over_Rp{i_alpha}, case_2.f{i_alpha} / case_2.Rx_over_Rb(i_alpha) / case_2.alpha(i_alpha) ./ abs( besselh(0,case_2.kRx{i_alpha})) ./ exp(-case_2.kh{i_alpha} * case_2.e2_over_h))
        hold on
    end
    legend(string(num2cell(case_2.alpha)))
    xlabel('kh R_x/R_p')
    ylabel('f/(\rho g \pi R_c^2 (R_x/R_b) \alpha H_0(k R_x) e^{-k e_2})')
    legend('location','best')

    figure
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
    figure
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

    figure
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

    figure
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

    figure
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

    figure
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

    figure
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

    figure
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
end

function fit_and_plot_case4(case_4, idxs_alpha_cell, idxs_beta_cell, ft, fo)
    % Fit and plot case 4 data
    x_string = 'k h R_x/R_p \alpha^2 / \beta';
    y_string = '$f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}\sqrt{kh})$';

    plot_fit_and_coeffs_sep_alpha_beta(case_4.kh_Rx_over_Rp_alpha2_over_beta, case_4.f_over_H0_exp_sqrt_kh, ...
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

    % version 1: semilogx with alpha^4 x-scaling
    fun = @(c, ia, ib) c * case_4.alpha(ia)^4;
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

    x_label_v1 = 'k h R_x/R_p \alpha^4 / \beta';
    y_label_v1 = '$\beta^2 \left[f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}kh R_x/R_b \alpha)\right]^2$ ';
    x_label_v2 = 'k h R_x/R_p \alpha^3 / \beta';
    y_label_v2 = '$\beta(1+\beta/\alpha) \left[f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}kh R_x/R_b \alpha)\right]^2$ ';

    % version 1: all alpha combined (manual + auto fits with coeffs)
    plot_fit_and_coeffs_all_alpha_beta(x_f1_v1, y_f1_v1, x_f3, y_base, y_pred_base, ...
        case_4.alpha, case_4.beta, @semilogx, x_label_v1, y_label_v1, [0 10], y_pred_f1_v1, ft, fo);

    % version 2: all alpha combined (manual + auto fits with coeffs)
    plot_fit_and_coeffs_all_alpha_beta(x_f1_v2, y_f1_v2, x_f3, y_base, y_pred_base, ...
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

    plot_fit_all_alpha_beta(x_new, y_new_x15, x_new, y_new_x15, y_pred_new, ...
        case_4.alpha, case_4.beta, @loglog, ...
        'k h R_x/R_p \alpha^3 / \beta', ...
        '$\beta(1+\beta/\alpha) \left[f/(\rho g \pi R_c^2  H_0(k R_x) e^{-k e_1}kh R_x/R_b \alpha)\right]^2 x^{1.5}$ ', ...
        [], y_pred_new);

    % Export fit data
    alpha_mat = cellfun(@(ia) case_4.alpha(ia), idxs_alpha_cell, 'UniformOutput', false);
    beta_mat  = cellfun(@(ib) case_4.beta(ib),  idxs_beta_cell,  'UniformOutput', false);
    case_4_xyab = [x_new(:) y_new(:) alpha_mat(:) beta_mat(:)];
    writecell(case_4_xyab, 'case_4_xyab.csv')
end

function plot_fit_and_coeffs_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string)
    % Plots manual fit, auto fit, and fit coefficients vs alpha/beta for separate-alpha subplots.
    plot_fit_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string, false);
    fits = plot_fit_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string, true);
    plot_fit_coeffs_vs_alpha_beta(fits, alpha_vec, beta_vec);
end

function plot_fit_and_coeffs_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
        alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell, ft, fo)
    % Plots manual fit, auto fit, and fit coefficients vs alpha/beta for all-alpha combined plots.
    plot_fit_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
        alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell);
    fits = plot_fit_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
        alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell, ft, fo);
    plot_fit_coeffs_vs_alpha_beta(fits, alpha_vec, beta_vec);
end

function plot_fit_coeffs_vs_alpha_beta(fits,alpha_vec,beta_vec)
    cols = {'r','g','b','k'};
    % coeffs vs beta for each alpha
    coeff_names = coeffnames(fits{1,1});
    figure
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
    figure
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

end

function fits = plot_fit_sep_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string, auto_fit_on)
    cols = {'r','g','b','k'};
    if auto_fit_on
        x_string = ['ln( ' x_string ')'];
        y_string = ['ln( ' y_string ')'];
    end
    fits = cell(size(alpha_vec));
    figure
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

function plot_fit_one_alpha_beta(x_cell, y_cell, alpha_vec, beta_vec, ft, fo, x_string, y_string)
    cols = {'r','g','b','k'};
    figure
    for i_alpha=1:3
        for j_beta = 1:4
            x = x_cell{i_alpha,j_beta};
            y = y_cell{i_alpha,j_beta};

            f = perform_auto_fit(x,y,ft,fo);
            h = plot(f,log(x),log(y));
            h(1).Color = cols{j_beta};
            h(2).Color = cols{j_beta};
            h(1).HandleVisibility = 'off';
            hold on
        end
    end
    legend(string(num2cell(beta_vec)),'location','best')
    xlabel(x_string)
    ylabel(y_string,'Interpreter','latex')
end

function fits = plot_fit_all_alpha_beta(x_f1_cell, y_f1_cell, x_f3_cell, y_cell, y_pred_cell, ...
        alpha_vec, beta_vec, f1_plot_fcn, x_f1_label, y_f1_label, y_f1_lim, y_f1_pred_cell, ft, fo)
    % Optional args:
    %   y_f1_pred_cell (12th): manual-fit predictions in y_f1 space for overlay on f1
    %   ft, fo         (13th, 14th): fittype/fitoptions for auto fitting; returns fits when provided
    has_f1_pred  = nargin >= 12 && ~isempty(y_f1_pred_cell);
    do_auto_fit  = nargin >= 14 && ~isempty(fo);
    fits = cell(length(alpha_vec), length(beta_vec));

    f1 = figure; f2 = figure; f3 = figure;
    for i_alpha = 1:length(alpha_vec)
        for j_beta = 1:length(beta_vec)
            alpha_beta_label = ['\alpha=' num2str(alpha_vec(i_alpha)) ', \beta=' num2str(beta_vec(j_beta))];
            x_f1 = x_f1_cell{i_alpha, j_beta};
            y_f1 = y_f1_cell{i_alpha, j_beta};

            figure(f1)
            if has_f1_pred
                % data as markers (no legend entry); manual fit as dashed (no legend entry)
                h = f1_plot_fcn(x_f1, y_f1, '*', 'HandleVisibility', 'off');
                hold on
                f1_plot_fcn(x_f1, y_f1_pred_cell{i_alpha, j_beta}, '--', ...
                    'Color', h.Color, 'HandleVisibility', 'off')
            else
                f1_plot_fcn(x_f1, y_f1, 'DisplayName', alpha_beta_label)
            end
            hold on

            if do_auto_fit
                f = perform_auto_fit(x_f1, y_f1, ft, fo);
                fits{i_alpha, j_beta} = f;
                % overlay auto fit as solid line with legend label
                x_fit_line = exp(linspace(log(min(x_f1(x_f1>0))), log(max(x_f1)), 100));
                y_fit_line = exp(feval(f, log(x_fit_line)));
                h_fit = f1_plot_fcn(x_fit_line, y_fit_line, 'DisplayName', alpha_beta_label);
                if has_f1_pred
                    h_fit.Color = h.Color;
                end
                y_pred_ij = exp(feval(f, log(x_f1)));
            else
                y_pred_ij = y_pred_cell{i_alpha, j_beta};
            end

            figure(f2)
            scatter(y_cell{i_alpha, j_beta}, y_pred_ij, 'DisplayName', alpha_beta_label)
            hold on

            figure(f3)
            scatter(x_f3_cell{i_alpha, j_beta}, y_pred_ij ./ y_cell{i_alpha, j_beta} - 1, 'DisplayName', alpha_beta_label)
            hold on
        end
    end

    figure(f1)
    legend('location','eastoutside')
    xlabel(x_f1_label)
    ylabel(y_f1_label,'Interpreter','latex')
    title('All case 4 data with fits')
    if ~isempty(y_f1_lim)
        ylim(y_f1_lim)
    end

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
