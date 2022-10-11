hydro = struct();
hydro = readWAMIT(hydro,'rm3.out',[]); % function from WECSim

w = hydro.w;
A = hydro.A(3,3,:);
A = A(:);
B = hydro.B(3,3,:);
B = B(:);
gamma = hydro.ex_ma(3,1,:);
gamma = gamma(:);

r = 10;
draft = 2;
g = 9.8;

k = w.^2 / g;

[A_MDOcean,B_MDOcean,gamma_MDOcean] = get_hydro_coeffs(r,k,draft);
A_MDOcean = ones(size(w))* A_MDOcean;

%% coeff comparison validation figure
figure
plot(w,A,'--',w,B,'--',w,gamma,'--','LineWidth',3,'HandleVisibility','off')
hold on
set(gca,"ColorOrderIndex",1)
plot(w,A_MDOcean,w,B_MDOcean,w,gamma_MDOcean)
xlim([0 1.5])
ylim([0 2500])
plot(0,0,'k-',0,0,'k--') % dummy plot so I can get 2 extra legend entries
legend('Added Mass A/\rho','Radiation Damping B/(\rho\omega)','Excitation Force \gamma/(\rhog)','Simulation (Analytical)','Actual (WAMIT BEM)')
title('Normalized Hydrodynamic Coefficients')
xlabel('Wave frequency \omega (rad/s)')
improvePlot

%% mean error
err_A = abs(A - A_MDOcean') ./ A;
err_B = abs(B - B_MDOcean') ./ B;
err_G = abs(gamma - gamma_MDOcean') ./ gamma;
mean(err_A(w<1.5))
mean(err_B(w<1.5))
mean(err_G(w<1.5))

%% R^2 table

R_A = corrcoef(A, A_MDOcean);
R_B = corrcoef(B, B_MDOcean);
R_G = corrcoef(gamma, gamma_MDOcean);

R2_A = R_A(1,2)^2
R2_B = R_B(1,2)^2
R2_G = R_G(1,2)^2

%% check B formula
figure
plot(w,B./gamma.^2*10000, w,w.^2/(2*9.8)*10000,'--')
legend('B/\gamma^2*\rho g^2/\omega', '\omega^2/2g')
xlim([0 1.5])
