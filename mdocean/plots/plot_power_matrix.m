function plot_power_matrix(X,p)

[~,~,P_matrix] = simulation(X, p);

P_matrix = P_matrix / 1e3; % convert W to kW
[CW_over_CW_max, P_wave, CW_max] = check_max_CW();%filename_uuid);

[T,Hs] = meshgrid(p.T,p.Hs);
figure

% subplot(1,3,1)
% contourf(T,Hs,P_matrix);
% xlabel('Wave Period T (s)')
% ylabel('Wave Height Hs (m)')
% title('Power (kW)')
% colorbar
% 
% subplot(1,3,2)
% contourf(T,Hs,p.JPD)
% xlabel('Wave Period T (s)')
% ylabel('Wave Height Hs (m)')
% title('Probability (%)')
% colorbar
% 
% subplot(1,3,3)
% contourf(T,Hs,P_matrix .* p.JPD/100)
% xlabel('Wave Period T (s)')
% ylabel('Wave Height Hs (m)')
% title('Weighted Power (kW)')
% colorbar

subplot(1,5,1)
contourf(T,Hs,P_wave);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Raw Wave Power Density (W/m)','FontSize',12)
colorbar

subplot(1,5,2)
contourf(T,Hs,CW_over_CW_max)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Device Capture Efficiency','FontSize',12)
colorbar

subplot(1,5,3)
contourf(T,Hs,CW_max)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Radiation Capture Width Limit (m)','FontSize',12)
colorbar

subplot(1,5,4)
contourf(T,Hs,p.JPD)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Probability at Site (%)','FontSize',12)
colorbar

subplot(1,5,5)
contourf(T,Hs,P_wave .* CW_over_CW_max .* CW_max .* p.JPD)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Weighted Power (W)','FontSize',12)
colorbar

%improvePlot

end
