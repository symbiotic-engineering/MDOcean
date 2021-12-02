function plot_power_matrix(X,p)

[~,~,~,~,~,~,~,~,~,P_matrix] = simulation(X, p);

[T,Hs] = meshgrid(p.T,p.Hs);
figure

subplot(1,3,1)
contourf(T,Hs,P_matrix);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Power (kW)')
colorbar

subplot(1,3,2)
contourf(T,Hs,p.JPD)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Probability (%)')
colorbar

subplot(1,3,3)
contourf(T,Hs,P_matrix .* p.JPD/100)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Weighted Power (kW)')
colorbar

improvePlot

end
