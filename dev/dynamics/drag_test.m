% script to show effect of drag and hydro coeffs
% compare against Yu and Li 2013 RANS paper, Figure 21
% see summer 24 research slideshow, slide 81

close all; clear 

reactive = true;

Cd = 0:1:2;

% singlebody
P_reactive = sweep_drag(true,false,Cd);
P_damping = sweep_drag(false,false,Cd);

% pct_diff = 100*(P_damping - P_reactive) ./ P_reactive

% multibody
sweep_drag(true,true,Cd)
sweep_drag(false,true,Cd)

function P_over_H2_2pt5 = sweep_drag(reactive,multibody,Cd)
    p = parameters();
    p.eff_pto = 1;
    p.use_multibody = multibody;
    
    if reactive
        p.control_type = 'Reactive';
        linestyle = '-';
    else
        p.control_type = 'Damping';
        linestyle = '--*';
    end

    if multibody
        bodies = 'Multibody';
    else
        bodies = 'Single Body';
    end

    b = var_bounds();
    X = [b.X_noms; 1];
    X(strcmp(b.var_names,'F_max')) = Inf;
    X(strcmp(b.var_names,'P_max')) = Inf;
    
    %p.T = 6.5;
    %p.Hs = 2.25;
    %p.JPD = 1;
    [~,H] = meshgrid(p.T, p.Hs);

    
    figure
    hold on
    for i = 1:length(Cd)
        p.C_d_float = Cd(i);
        disp(Cd(i))
        [~, P_matrix, ~, val] = simulation(X,p);
        P_over_H2 = P_matrix ./ H.^2;
        P_over_H2_2pt5 = P_over_H2(H==2.25);
        %P_over_H2_2pt5 = .5 * (P_over_H2(H==2.25) + P_over_H2(H==2.75));
        plot(p.T, P_over_H2_2pt5/1000,linestyle,'DisplayName',['Cd=' num2str(Cd(i))])

        disp(val.X_f(H==2.25))
        disp(val.force_ptrain)
    end
    legend
    xlabel('T (s)')
    ylabel('P / H^2 (kW/m^2) at H=2.5m')
    
    T_new = [4 4.5 p.T]; % extend to lower periods
    w = 2*pi./T_new;
    k = w.^2 / p.g;
    lambda = 2*pi./k;
    P_over_H2_max = p.rho_w * p.g^2 * T_new .* lambda / (64 * pi^2);
    
    plot(T_new,P_over_H2_max/1000,...
        'DisplayName','Theoretical limit')
    ylim([0 185])
    xlim([4 18])
    title([p.control_type ' Control, ' bodies])
    improvePlot
end