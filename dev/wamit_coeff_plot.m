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
fudge_factor =  5;
draft = 2;
g = 9.8;
rho_w = 1000;

k = w.^2 / g;
V_g = g ./(2*w);

A_MDOcean = ones(size(w))* 1/2 * 4/3 * pi * r^3 * 0.63;
r_k_term = r^2 - 1/8 * k.^2 * r^4 + 1/192 * k.^4 * r^6 - 1/9216 * k.^6 * r^8;
gamma_MDOcean   = pi * exp(-k * draft * fudge_factor) .* r_k_term; 
B_MDOcean = k/2 .* gamma_MDOcean.^2;
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
