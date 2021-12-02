
x = linspace(0,2*pi,100);
y = 1.5*sin(x);
idx_sat = abs(y)>1;
ysat = y;
ysat(idx_sat) = 1*sign(y(idx_sat));

yfund = 1.18*sin(x);

figure
h = plot(x,y,'r--');
hold on
plot(x,ysat,'r',x,yfund,'b')
improvePlot
set(h,'HandleVisibility','off')

legend('Saturated Signal','Fundamental Amplitude')

