close all
years = 10;

capex = 10;
opex = 1;
costs = [-capex repmat(-opex,1,years)];

FCR = 0.113;
levelized_cost = FCR * capex + opex;

revenue = 2;
revenues = [0 repmat(revenue,1,years)];

yrange = 1.1*[-capex,revenue];

% just cost
figure('Color','w')
subplot 121
bar(costs,'r')
axis off
ylim(yrange)
xlims = xlim;

subplot 122
bar(-levelized_cost,'r')
axis off
ylim(yrange)
xlim(xlims)

% cost and revenue
figure('Color','w')
subplot 121
bar(costs,'r')
hold on
bar(revenues,'g')
axis off
ylim(yrange)

subplot 122
bar(-levelized_cost,'r')
hold on
bar(revenue,'g')
axis off
ylim(yrange)
xlim(xlims)


