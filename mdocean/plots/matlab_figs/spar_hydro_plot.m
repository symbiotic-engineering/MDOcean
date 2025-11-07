function f = spar_hydro_plot()

    r = linspace(0, 1);
    root_r = sqrt(1 - r.^2);
    A_s_over_rho_Dd_3 = 1 / 3 - 1 / 4 * r.^2 .* root_r - 1 / 12 * (1 - root_r).^2 .* (2 + root_r);

    f = figure;
    plot(r, A_s_over_rho_Dd_3 * 3);
    xlabel('$r=\frac{D_s}{D_d}$', 'Interpreter', 'latex');
    ylabel('$A_s/(\frac{1}{3}\rho_w D_d^3)$', 'Interpreter', 'latex');
    title('Spar Added Mass');
    improvePlot;
    f.Position(3:4) = [550 380];

end
