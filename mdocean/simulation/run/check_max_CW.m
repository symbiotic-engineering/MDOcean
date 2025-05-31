function [hydro_ratio, P_wave, CW_max, ...
    P_elec, force_sat_ratio, drag_ratio, eff] = check_max_CW(filename_uuid, X, p, b, plot_on)
    
    if nargin<4
        p = parameters();
        b = var_bounds();
        if nargin==0
            filename_uuid = '';
        end
        b.filename_uuid = filename_uuid;
    
        [X,val,~,P_elec] = max_avg_power(p,b);    
    else
        [~,~,P_elec,~,val] = simulation(X, p);
    end
    if nargin<5
        plot_on = true;
    end
    P_mech = val.P_mech;
    force_sat_ratio = val.P_sat_ratio;
    P_mech_unsat =  P_mech ./ force_sat_ratio;
    eff = P_elec ./ P_mech;

    % solution with no force saturation
    X_unsat = X;
    idx_F = strcmp(b.var_names,'F_max');
    X_unsat(idx_F) = Inf;

    % solution with no force sat and no drag
    p_no_drag = p;
    p_no_drag.C_d_float = 0;
    p_no_drag.C_d_spar  = 0;
    [~,~,~,~,val_no_drag] = simulation(X_unsat, p_no_drag);
    P_no_drag = val_no_drag.P_mech;
    drag_ratio = P_mech_unsat ./ P_no_drag;

    % calculate capture width of unsat
    [T,Hs] = meshgrid(p.T,p.Hs);
    P_wave = p.rho_w * p.g^2 / (64*pi) * T .* Hs.^2;
    CW_unsat = P_no_drag ./ P_wave;

    % compare to maximum capture width
    CW_max = p.g * T.^2 / (4*pi^2);
    hydro_ratio = CW_unsat ./ CW_max;


    if plot_on
        figure
        contourf(T,Hs,hydro_ratio)
        xlabel('Wave Period T (s)')
        ylabel('Wave Height Hs (m)')
        colorbar
        grid on
        title('CW / CW_{max}')
    
        figure
        plot(p.T,hydro_ratio(1,:))
        xlabel('Wave Period T (s)')
        ylabel('CW / CW_{max}')
        hold on
        plot([min(p.T) max(p.T)], [1 1], 'k--')
        improvePlot
    end
end

