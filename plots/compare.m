close all

nominal	= [20	0.3	.2	.8 5	10	0.79 1];
cheapest = [13.7	0.44	.1 .79 4.5 14 .85 1];
minVariation = [12 .5 .1 .98 11 160 36.6 1];
balanced = [13.0468, 0.4599, 0.1000, 0.9738, 12.7445, 30.6899, 10.1115, 1.0000];


p = parameters();

X = [nominal;cheapest;minVariation;balanced];

col = {'k','b','r','g'};

%% geometry comparison
figure
for i=1:4
    x = X(i,:);
    hold on
    visualize_geometry(x,p,false,col{i})
end
legend('Nominal','Min LCOE','Min Variation','Balanced')

%% power probability comparison
figure
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
ylabel(t,'Probability (-)','FontWeight','bold')
for i=1:4
    x = X(i,:);
    hold on
    t = power_PDF(x,p,col{i},t);
end
legend('Nominal','Min LCOE','Min Variation','Balanced','location','best')

%% power matrix comparison
titles = {'Nominal','Min LCOE','Min Variation','Balanced'};
figure
for i=1:4
    x = X(i,:);   
    [~,~,~,~,~,~,~,~,~,P_matrix] = simulation(x, p);
    P_matrix = P_matrix / 1e3; % convert W to kW
    [T,Hs] = meshgrid(p.T,p.Hs); 
    
    subplot(2,2,i);
    hold on
    contourf(T,Hs,P_matrix .* p.JPD/100,[0 5 10:10:60])
    improvePlot
    grid on
    title(titles{i})
    caxis([0 60])
    ax = gca;
    if i==2 || i==4
        xoffset = -ax.Position(3)*.25;
    else
        xoffset = 0;
    end
    if i==3 || i==4
        yoffset = ax.Position(4)*.2;
    else
        yoffset = 0;
    end
    set(ax,'Position',[ax.Position(1)+xoffset ax.Position(2)+yoffset ...
        ax.Position(3)*.8 ax.Position(4)*.8])
    set(ax.Title,'FontWeight','normal')
    set(ax,'Layer','top','GridColor',[.9 .9 .9])
end
cb = colorbar('Position',[.8 .175 .08 .685]);
%cb.Label.String = 'Weighted Power (kW)';
sgtitle('Weighted Power (kW)','FontWeight','bold','FontSize',20)
text(-10,-2,'Wave Period T (s)','FontWeight','bold','FontSize',16)
text(-30,5,'Wave Height Hs (m)','FontWeight','bold','FontSize',16,'Rotation',90)
