function intermed_result_struct = analysis_fcn(p, b)
% analysis_fcn  Run the geometry sweep simulation loop.
%
% Sweeps six geometric ratios (a_1/a_2, a_2/h, d_1/h, d_2/d_1, a_3/a_1)
% and water depth h, computing the hydro ratio CW/CW_max and the wave
% parameter m_0*h at every grid point.
%
% :param p: Parameter struct (from parameters())
% :param b: Design variable bounds struct (from var_bounds())
% :returns: intermed_result_struct  Struct with hydro_ratio_result,
%           m0h_stored, val (simulation outputs), all grid arrays, the
%           sweep parameter vectors, and p.

    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')

    n = 3; % number of points per dimension

    zero_to_one   = linspace(0.01, 0.99, n);
    zero_to_three = linspace(0.01, 3,    n);

    h_vec = [10 200]; % gives a range of m0h of 0.36-2.06 and 2.40-39.79 for default wave periods
    a1_a2 = zero_to_one;
    a2_h  = zero_to_three;
    d1_h  = zero_to_one;
    d2_d1 = zero_to_one;
    a3_a1 = 1 + zero_to_three;

    [H, A1_A2, A2_H, D1_H, D2_D1, A3_A1] = ndgrid(h_vec, a1_a2, a2_h, d1_h, d2_d1, a3_a1);

    A1_H = A1_A2 .* A2_H;
    D2_H = D2_D1 .* D1_H;
    A3_H = A3_A1 .* A1_H;

    A1 = A1_H .* H;
    A2 = A2_H .* H;
    A3 = A3_H .* H;
    D1 = D1_H .* H;
    D2 = D2_H .* H;

    X = [b.X_noms; 1];

    m0h_stored         = zeros([length(p.T), size(A1)]);
    hydro_ratio_result = zeros([length(p.T), size(A1)]);

    for i = 1:numel(A1)
        disp(['Running ' num2str(i) ' of ' num2str(numel(A1))])

        D_s   = 2 * A1(i);
        D_f   = 2 * A2(i);
        D_d   = 2 * A3(i);
        T_s   = D1(i);
        T_f_2 = D2(i);

        % Construct per-iteration parameter struct and design vector
        p_i = p;
        p_i.h              = H(i);
        p_i.T_s_over_D_s   = T_s   / D_s;
        p_i.D_d_over_D_s   = D_d   / D_s;
        X_i = X;
        X_i(strcmp(b.var_names, 'D_s'))   = D_s;
        X_i(strcmp(b.var_names, 'D_f'))   = D_f;
        X_i(strcmp(b.var_names, 'h_s'))   = T_s + 5;
        X_i(strcmp(b.var_names, 'T_f_2')) = T_f_2;

        subs = {ind2sub(size(A1), i)};
        m0h_stored(:, subs{:}) = dispersion(2*pi ./ p.T, p_i.h, p.g) * p_i.h;

        % Run simulation (drag and force/power sat are turned off inside check_max_CW)
        [hydro_ratio, ~, ~, ~, ~, ~, ~, out] = check_max_CW('', p_i, X_i, false);
        if i == 1
            val = out;
        else
            val(i) = out;
        end

        hydro_ratio_result(:, subs{:}) = max(hydro_ratio);
    end

    % --- Save intermediate results ---
    % Main simulation outputs
    intermed_result_struct.hydro_ratio_result = hydro_ratio_result;
    intermed_result_struct.m0h_stored         = m0h_stored;
    intermed_result_struct.val                = val;

    % Dimensional grid arrays (needed for CWR computation)
    intermed_result_struct.H     = H;
    intermed_result_struct.A1_A2 = A1_A2;
    intermed_result_struct.A2_H  = A2_H;
    intermed_result_struct.D1_H  = D1_H;
    intermed_result_struct.D2_D1 = D2_D1;
    intermed_result_struct.A3_A1 = A3_A1;
    intermed_result_struct.A1    = A1;
    intermed_result_struct.A2    = A2;
    intermed_result_struct.D1    = D1;
    intermed_result_struct.D2    = D2;

    % Sweep parameter vectors (needed for plot annotations)
    intermed_result_struct.h_vec  = h_vec;
    intermed_result_struct.a1_a2  = a1_a2;
    intermed_result_struct.a2_h   = a2_h;
    intermed_result_struct.d1_h   = d1_h;
    intermed_result_struct.d2_d1  = d2_d1;
    intermed_result_struct.a3_a1  = a3_a1;

    % Parameter struct (p.T and p.g needed for CWR in post-processing)
    intermed_result_struct.p = p;
end
