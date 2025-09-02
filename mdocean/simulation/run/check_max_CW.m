function [hydro_ratio, P_wave, CW_max, ...
    P_elec, force_sat_ratio, drag_ratio, eff] = check_max_CW(filename_uuid, p, b_or_X, plot_on)
% Checks whether the maximum radiation-limited capture width is violated.
% b_or_X: provide b (or no input) if you want to use the design with highest avg power (worst case),
% vs provide X if you want to check a certain design.
    
    if nargin<1
        filename_uuid = '';
    end
    if nargin<2
        p = parameters();
    end
    if nargin<3
        b_or_X = var_bounds();
    end
    if nargin<4
        plot_on = true;
    end
    
    if isstruct(b_or_X)
        b = b_or_X;
        b.filename_uuid = filename_uuid;
    
        [X,val,~,P_elec] = max_avg_power(p,b);    
    elseif isfloat(b_or_X)
        X = b_or_X;
        [~,P_elec,~,val] = simulation(X, p);
    else
        msg = ['b_or_X is ' class(b_or_X) ...
            ' but should be either struct (b) or float (X)'];
        error(msg)
    end

    P_mech = val.P_mech;
    force_sat_ratio = val.P_sat_ratio;
    P_mech_unsat =  P_mech ./ force_sat_ratio;
    eff = P_elec ./ P_mech;

    % solution with no force sat, power sat, or drag
    p_no_drag_no_sat = p;
    p_no_drag_no_sat.C_d_float = 0;
    p_no_drag_no_sat.C_d_spar  = 0;
    p_no_drag_no_sat.use_force_sat = 0;
    p_no_drag_no_sat.use_power_sat = 0;
    [~,~,~,val_no_drag] = simulation(X, p_no_drag_no_sat);
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

