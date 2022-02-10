close all

nominal	= [20	0.3	30	1 10	10	0.8];
cheapest = [7	0.4	12	1 100	1	0.3];
minVariation = [	21	0.5	32	1 81	11	1.2];
balanced = [	6	1	33	1 62	8	0.4];


p = parameters();

X = [nominal;cheapest;minVariation;balanced];

col = {'k','b','r','g'};

figure
for i=1:4
    x = X(i,:);
    hold on
    visualize_geometry(x,p,false,col{i})
end
legend('Nominal','Min LCOE','Min Variation','Balanced')

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
