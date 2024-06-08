function [mean_abs_err, R2] = hydro_coeff_err(plot_on)

if nargin==0
    plot_on = true;
end

hydro = struct();
hydro = readWAMIT(hydro,'rm3.out',[]); % function from WECSim

w_max = 1.5;
w = hydro.w;
w = w(w<w_max);
A = hydro.A(3,3,w<w_max);
A = A(:);
B = hydro.B(3,3,w<w_max);
B = B(:);
gamma = hydro.ex_ma(3,1,w<w_max);
gamma = gamma(:);

a2 = 10;
a1 = 3;
d2 = 2;
d1 = 35;
g = 9.8;
h = 100;

k = w.^2 / g;

use_MEEM = true;
if use_MEEM
    [A_MDOcean,B_MDOcean,gamma_MDOcean] = get_hydro_coeffs_MEEM(a2, k, d2, a1, d1, h); 
else
    [A_MDOcean,B_MDOcean,gamma_MDOcean] = get_hydro_coeffs(a2, k, d2);
    A_MDOcean = ones(size(w))* A_MDOcean;
end


%% mean error
err_A = abs(A - A_MDOcean') ./ A;
err_B = abs(B - B_MDOcean') ./ B;
err_G = abs(gamma - gamma_MDOcean') ./ gamma;
mA = mean(err_A);
mB = mean(err_B);
mG = mean(err_G);

mean_abs_err = [mA mB mG];

%% R^2
R_A = corrcoef(A, A_MDOcean);
R_B = corrcoef(B, B_MDOcean);
R_G = corrcoef(gamma, gamma_MDOcean);

R2_A = R_A(1,2)^2;
R2_B = R_B(1,2)^2;
R2_G = R_G(1,2)^2;

R2 = [R2_A R2_B R2_G];

if plot_on
    %% coeff comparison validation figure
    figure
    plot(w,A,'--',w,B,'--',w,gamma,'--','LineWidth',3,'HandleVisibility','off')
    hold on
    set(gca,"ColorOrderIndex",1)
    plot(w,A_MDOcean,w,B_MDOcean,w,gamma_MDOcean)
    ylim([0 3500])
    plot(0,0,'k-',0,0,'k--') % dummy plot so I can get 2 extra legend entries
    legend('Added Mass A/\rho','Radiation Damping B/(\rho\omega)','Excitation Force \gamma/(\rhog)','Simulation (Analytical)','Actual (WAMIT BEM)')
    title('Normalized Hydrodynamic Coefficients')
    xlabel('Wave frequency \omega (rad/s)')
    improvePlot
    
    %% check B formula
    figure
    plot(w,B./gamma.^2*10000, w,w.^2/(2*9.8)*10000,'--')
    legend('B/\gamma^2*\rho g^2/\omega', '\omega^2/2g')
end
