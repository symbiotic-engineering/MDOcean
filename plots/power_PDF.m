function t = power_PDF(X,p,color,t)

[~,~,~,~,~,~,~,~,~,P_matrix] = simulation(X, p);
P_matrix = P_matrix / 1e3; % convert W to kW
P_counts = repelem(P_matrix(:), round(p.JPD(:)*100));

% PDF
ax1 = nexttile(1);
%h = histogram(P_counts,'BinWidth',50,'Normalization','pdf');
[N,edges] = histcounts(P_counts,'BinWidth',100,'Normalization','pdf');
edge_centers = 1/2 *(edges(1:end-1) + edges(2:end));
h = plot(edge_centers,N);
grid on
if nargin > 2
    h.Color = color;
end
hold on
title('Probability Density Function')
yy = ylabel('Probability (-)');
yy.Position = [-100 -.001 0];


% CDF
ax2 = nexttile(2);
h = cdfplot(P_counts);
hold on
if nargin > 2
    h.Color = color;
end
title('Cumulative Distribution Function')
ylabel('')

linkaxes([ax1,ax2],'x');
xlabel('Power (kW)')
xticklabels(ax1,{})
xticks(ax1,xticks(ax2))
xlim([0 1000])
improvePlot

nexttile(1)


end