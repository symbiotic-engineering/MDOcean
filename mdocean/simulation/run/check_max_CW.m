function [ratio, P_wave, CW_max] = check_max_CW()%filename_uuid)
    p = parameters();
    b = var_bounds();
    %b.filename_uuid = filename_uuid;

    [~,val] = max_avg_power(p,b);

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

