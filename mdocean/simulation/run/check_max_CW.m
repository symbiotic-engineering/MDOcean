function [ratio] = check_max_CW()
    p = parameters();
    p.cost_m = [0 0 0]; % hack that makes cost constant, so minimizing LCOE is actually maximizing power
    p.power_max = Inf;
    
    b = var_bounds();
    x0 = b.X_start_struct;

    % run LCOE minimization (effectively power maximization due to hack above)
    X_opt = gradient_optim(x0,p,b,1); 

    % plug back into simulation to get unsaturated power
    [~, ~, ~, ~, val] = simulation(X_opt,p);
    P_mech = val.P_mech;

    % calculate capture width
    [T,Hs] = meshgrid(p.T,p.Hs);
    P_wave = p.rho_w * p.g^2 / (64*pi) * T .* Hs.^2;
    CW = P_mech ./ P_wave;

    % compare to maximum capture width
    CW_max = p.g * T.^2 / (4*pi^2);

    ratio = CW ./ CW_max;

    figure
    contourf(T,Hs,ratio)
    xlabel('Wave Period T (s)')
    ylabel('Wave Height Hs (m)')
    colorbar
    grid on
    title('CW / CW_{max}')

    figure
    plot(p.T,ratio(1,:))
    xlabel('Wave Period T (s)')
    ylabel('CW / CW_{max}')
    hold on
    plot([min(p.T) max(p.T)], [1 1], 'k--')
    improvePlot
end

