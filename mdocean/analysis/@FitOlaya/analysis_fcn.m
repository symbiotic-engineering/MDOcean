function intermed_result_struct = analysis_fcn(~,~)
% Function analysis_fcn
%
% :param ~: ~
% :param ~: ~
% :returns: Intermediate results struct with raw data and derived signals

    raw = load_raw_data();
    [case_2, case_4, idxs_alpha_cell, idxs_beta_cell] = create_signals(raw);

    intermed_result_struct.raw             = raw;
    intermed_result_struct.case_2          = case_2;
    intermed_result_struct.case_4          = case_4;
    intermed_result_struct.idxs_alpha_cell = idxs_alpha_cell;
    intermed_result_struct.idxs_beta_cell  = idxs_beta_cell;
end

function raw = load_raw_data()
    % Load raw digitized data from CSV files for figures 3b, 4b, 12a-12d
    this_dir = fileparts(mfilename('fullpath'));
    data_dir = fullfile(this_dir, '..', '..', '..', 'dev', 'hydro_coeffs', 'damping_plate');

    data_3b_alpha_2    = readmatrix(fullfile(data_dir, '3b_red.csv'));
    data_3b_alpha_1    = readmatrix(fullfile(data_dir, '3b_black.csv'));
    data_3b_alpha_pt51 = readmatrix(fullfile(data_dir, '3b_blue.csv'));

    raw.kRb_3b_alpha_2    = data_3b_alpha_2(:,1);    raw.f_3b_alpha_2    = data_3b_alpha_2(:,2);
    raw.kRb_3b_alpha_1    = data_3b_alpha_1(:,1);    raw.f_3b_alpha_1    = data_3b_alpha_1(:,2);
    raw.kRb_3b_alpha_pt51 = data_3b_alpha_pt51(:,1); raw.f_3b_alpha_pt51 = data_3b_alpha_pt51(:,2);

    data_4b_alpha_8 = readmatrix(fullfile(data_dir, '4b_red.csv'));
    data_4b_alpha_2 = readmatrix(fullfile(data_dir, '4b_blue.csv'));
    data_4b_alpha_1 = readmatrix(fullfile(data_dir, '4b_black.csv'));

    raw.kRb_4b_alpha_8 = data_4b_alpha_8(:,1); raw.f_4b_alpha_8 = data_4b_alpha_8(:,2);
    raw.kRb_4b_alpha_2 = data_4b_alpha_2(:,1); raw.f_4b_alpha_2 = data_4b_alpha_2(:,2);
    raw.kRb_4b_alpha_1 = data_4b_alpha_1(:,1); raw.f_4b_alpha_1 = data_4b_alpha_1(:,2);

    data_12a_alpha_1    = readmatrix(fullfile(data_dir, '12a_red.csv'));
    data_12a_alpha_pt75 = readmatrix(fullfile(data_dir, '12a_blue.csv'));
    data_12a_alpha_pt51 = readmatrix(fullfile(data_dir, '12a_black.csv'));

    raw.kRb_12a_alpha_1    = data_12a_alpha_1(:,1);    raw.f_12a_alpha_1    = data_12a_alpha_1(:,2);
    raw.kRb_12a_alpha_pt75 = data_12a_alpha_pt75(:,1); raw.f_12a_alpha_pt75 = data_12a_alpha_pt75(:,2);
    raw.kRb_12a_alpha_pt51 = data_12a_alpha_pt51(:,1); raw.f_12a_alpha_pt51 = data_12a_alpha_pt51(:,2);

    data_12b_alpha_1    = readmatrix(fullfile(data_dir, '12b_red.csv'));
    data_12b_alpha_pt75 = readmatrix(fullfile(data_dir, '12b_blue.csv'));
    data_12b_alpha_pt51 = readmatrix(fullfile(data_dir, '12b_black.csv'));

    raw.kRb_12b_alpha_1    = data_12b_alpha_1(:,1);    raw.f_12b_alpha_1    = data_12b_alpha_1(:,2);
    raw.kRb_12b_alpha_pt75 = data_12b_alpha_pt75(:,1); raw.f_12b_alpha_pt75 = data_12b_alpha_pt75(:,2);
    raw.kRb_12b_alpha_pt51 = data_12b_alpha_pt51(:,1); raw.f_12b_alpha_pt51 = data_12b_alpha_pt51(:,2);
    raw.f_12b_alpha_pt75(raw.f_12b_alpha_pt75<0) = 1e-3;

    data_12c_alpha_1    = readmatrix(fullfile(data_dir, '12c_red.csv'));
    data_12c_alpha_pt75 = readmatrix(fullfile(data_dir, '12c_blue.csv'));
    data_12c_alpha_pt51 = readmatrix(fullfile(data_dir, '12c_black.csv'));

    raw.kRb_12c_alpha_1    = data_12c_alpha_1(:,1);    raw.f_12c_alpha_1    = data_12c_alpha_1(:,2);
    raw.kRb_12c_alpha_pt75 = data_12c_alpha_pt75(:,1); raw.f_12c_alpha_pt75 = data_12c_alpha_pt75(:,2);
    raw.kRb_12c_alpha_pt51 = data_12c_alpha_pt51(:,1); raw.f_12c_alpha_pt51 = data_12c_alpha_pt51(:,2);

    data_12d_alpha_1    = readmatrix(fullfile(data_dir, '12d_red.csv'));
    data_12d_alpha_pt75 = readmatrix(fullfile(data_dir, '12d_blue.csv'));
    data_12d_alpha_pt51 = readmatrix(fullfile(data_dir, '12d_black.csv'));

    raw.kRb_12d_alpha_1    = data_12d_alpha_1(:,1);    raw.f_12d_alpha_1    = data_12d_alpha_1(:,2);
    raw.kRb_12d_alpha_pt75 = data_12d_alpha_pt75(:,1); raw.f_12d_alpha_pt75 = data_12d_alpha_pt75(:,2);
    raw.kRb_12d_alpha_pt51 = data_12d_alpha_pt51(:,1); raw.f_12d_alpha_pt51 = data_12d_alpha_pt51(:,2);
    raw.f_12d_alpha_pt51(raw.f_12d_alpha_pt51<0) = 1e-3;
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
