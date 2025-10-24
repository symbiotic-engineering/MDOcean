function figs = sin_desc_fcn_demo()

%% force saturation
x = linspace(0,1,100);
y = 1.5*sin(x*2*pi);
idx_sat = abs(y)>1;
ysat = y;
ysat(idx_sat) = 1*sign(y(idx_sat));

yfund = 1.18*sin(x*2*pi);

h1 = figure;
h = plot(x,y,'r--');
hold on
plot(x,ysat,'r',x,yfund,'b')
xlabel('Normalized Time (-)')
ylabel('Normalized Powertrain Force (-)')
improvePlot
set(h,'HandleVisibility','off')

legend('Saturated Signal','Fundamental Amplitude','Interpreter','latex')

% see what happens when added to phase shifted signal
shift_on = false;
if shift_on
    shift = 33;
    y_shift = .5/1.5 * [y(shift:end) y(1:(shift-1))];
    plot(x,y_shift,x,ysat+y_shift,x,yfund+y_shift)
end

% alpha plot

h2 = figure;
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

%% drag
drag_DF = 8/(3*pi)*sin(x*2*pi);
drag_actual = sin(x*2*pi) .* abs(sin(x*2*pi));

h3 = figure;
plot(x,drag_actual,'r',x,drag_DF,'b')
xlabel('Normalized Time (-)')
ylabel('Normalized Drag Force (-)')
legend('Squared Signal $\sin(\omega t)|\sin(\omega t)|$',...
    'Fundamental Amplitude $\frac{8}{3\pi}\sin(\omega t)$','Interpreter','latex')
improvePlot

% Return handles in creation order
figs = [h1, h2, h3];

end