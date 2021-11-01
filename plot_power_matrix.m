function plot_power_matrix(X,p)

[~,~,~,~,~,~,~,~,P_matrix,~] = simulation(X, p);

[T,Hs] = meshgrid(p.T,p.Hs);
figure

subplot(1,3,1)
contourf(T,Hs,P_matrix);
xlabel('T')
ylabel('Hs')
title('Power Matrix')
colorbar

subplot(1,3,2)
contourf(T,Hs,p.JPD)
xlabel('T')
ylabel('Hs')
title('Joint Probability Distribution')
colorbar

subplot(1,3,3)
contourf(T,Hs,P_matrix .* p.JPD/100)
xlabel('T')
ylabel('Hs')
title('Probability Weighted Power Matrix')
colorbar

end
