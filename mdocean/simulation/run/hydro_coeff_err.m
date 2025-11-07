function [mean_abs_err, R2, fig] = hydro_coeff_err(plot_on)

    if nargin == 0
        plot_on = true;
    end

    hydro = struct();
    hydro = readWAMIT(hydro, 'rm3.out', []); % function from WECSim

    w_max = 1.5;
    w = hydro.w;
    w = w(w < w_max);
    A = hydro.A(3, 3, w < w_max);
    A = A(:);
    B = hydro.B(3, 3, w < w_max);
    B = B(:);
    gamma = hydro.ex_ma(3, 1, w < w_max);
    gamma = gamma(:);
    gamma_ph = -hydro.ex_ph(3, 1, w < w_max);
    gamma_ph = gamma_ph(:);

    a2 = 10;
    a1 = 3;
    d2 = 2;
    d1 = 35;
    g = 9.8;
    h = 300;
    a3 = 15;

    spar_excitation_coeffs = get_spar_exc(g);

    k = w.^2 / g;
    harmonics = 20;

    use_MEEM = true;
    if use_MEEM
        [A_MDOcean, ~, ~, B_MDOcean, ~, ~, ...
         gamma_MDOcean, ~, gamma_ph_MDOcean, ~] = get_hydro_coeffs_MEEM(a2, k, d2, a1, d1, a3, h, g, w, harmonics, spar_excitation_coeffs);
    else
        [A_MDOcean, B_MDOcean, gamma_MDOcean] = get_hydro_coeffs(a2, k, d2);
        A_MDOcean = ones(size(w)) * A_MDOcean;
    end

    %% mean error
    err_A = abs(A - A_MDOcean') ./ A;
    err_B = abs(B - B_MDOcean') ./ B;
    err_G = abs(gamma - gamma_MDOcean') ./ gamma;
    err_Gph = abs(gamma_ph - gamma_ph_MDOcean') ./ gamma_ph;
    err_Gph(isnan(err_Gph)) = 0;

    mA = mean(err_A);
    mB = mean(err_B);
    mG = mean(err_G);
    mGph = mean(err_Gph);

    mean_abs_err = [mA mB mG mGph];

    %% R^2
    R_A = corrcoef(A, A_MDOcean);
    R_B = corrcoef(B, B_MDOcean);
    R_G = corrcoef(gamma, gamma_MDOcean);
    R_Gph = corrcoef(gamma_ph, gamma_ph_MDOcean);

    R2_A = R_A(1, 2)^2;
    R2_B = R_B(1, 2)^2;
    R2_G = R_G(1, 2)^2;
    R2_Gph = R_Gph(1, 2)^2;

    R2 = [R2_A R2_B R2_G R2_Gph];

    if plot_on
        %% coeff comparison validation figure
        fig = figure;
        hold on;
        colororder({'k', 'k'});

        blue = [0 0.4470 0.7410];
        red = [0.8500 0.3250 0.0980];
        yellow = [0.9290 0.6940 0.1250];
        purple = [0.4940 0.1840 0.5560];

        yyaxis left;
        plot(w, A, '--', 'color', blue, 'Linewidth', 3, 'HandleVisibility', 'off');
        plot(w, B, '--', 'color', red, 'LineWidth', 3, 'HandleVisibility', 'off');
        plot(w, gamma, '--', 'color', yellow, 'LineWidth', 3, 'HandleVisibility', 'off');
        yyaxis right;
        plot(w, gamma_ph, '--', 'color', purple, 'LineWidth', 3, 'HandleVisibility', 'off');
        set(gca, "ColorOrderIndex", 1);

        yyaxis left;
        ylim([0, 4000]);
        plot(w, A_MDOcean, '-', 'color', blue, 'LineWidth', 3);
        plot(w, B_MDOcean, '-', 'color', red, 'Linewidth', 3);
        plot(w, gamma_MDOcean, '-', 'color', yellow, 'Linewidth', 3);
        ylabel('A/\rho (m^{3}), B/(\rho\omega) (m^{3}/(rad^{2}s^{2})), |\gamma|/(\rhog) (m^{2})');
        yyaxis right;
        ylim([-1.8, 1.1]);
        ylabel('\angle\gamma (rad)');
        plot(w, gamma_ph_MDOcean, '-', 'color', purple, 'Linewidth', 3);

        plot(0, 0, 'k-', 0, 0, 'k--'); % dummy plot so I can get 2 extra legend entries
        leg = legend('Added Mass A/\rho', 'Radiation Damping B/(\rho\omega)', 'Excitation Force Magnitude |\gamma|/(\rhog)', ...
                     'Excitation Phase \angle\gamma', ...
                     'MDOcean Simulation (Semi-Analytical MEEM)', '"Ground Truth" Simulation (WAMIT BEM)');
        title('Normalized Hydrodynamic Coefficients');
        xlabel('Wave frequency \omega (rad/s)');
        improvePlot;
        set(fig, 'Position', [100 100 697.8 600]);
        set(leg, 'Position', [0.1692 0.6372 0.7027 0.2615]);

        %% check B formula
        %     figure
        %     plot(w,B./gamma.^2*10000, w,w.^2/(2*9.8)*10000,'--')
        %     legend('B/\gamma^2*\rho g^2/\omega', '\omega^2/2g')
    else
        fig = gobjects(1);
    end
