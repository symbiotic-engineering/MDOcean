%close all

p = parameters();
b = var_bounds(p);
x0_input = b.X_start_struct;

Xs_opt = gradient_optim(x0_input,p,b);
cheapest = Xs_opt(:,1);
minVariation = Xs_opt(:,2);

p.LCOE_max = 0.2;
Xs_opt = gradient_optim(x0_input,p,b);
balanced = Xs_opt(:,2);

nominal	= [b.X_noms' 1];
%%
X = [nominal;cheapest';minVariation';balanced'];

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
    [~,~,~,~,~,~,~,~,~,~,P_matrix] = simulation(x, p);
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
