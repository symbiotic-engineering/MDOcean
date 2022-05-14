
x = linspace(0,1,100);
y = 1.5*sin(x*2*pi);
idx_sat = abs(y)>1;
ysat = y;
ysat(idx_sat) = 1*sign(y(idx_sat));

yfund = 1.18*sin(x*2*pi);

figure
h = plot(x,y,'r--');
hold on
plot(x,ysat,'r',x,yfund,'b')
xlabel('Normalized Time')
ylabel('Normalized Powertrain Force')
improvePlot
set(h,'HandleVisibility','off')

legend('Saturated Signal','Fundamental Amplitude')

figure
r = linspace(1,8,50); % r = Fp / Fmax
alpha = 2/pi * ( r.*asin(1./r) + sqrt(1 - 1./r.^2) );
plot(r,alpha)
hold on
plot(r,4/pi*ones(size(r)),'k--')
text(1.25,1.26,'4/\pi','FontSize',18)
xlim([1 8])
xticks(1:8)
yticks(1:.1:1.3)
xlabel('Force Ratio F_p/F_{max}')
ylabel('Fundamental Ratio \alpha')
improvePlot