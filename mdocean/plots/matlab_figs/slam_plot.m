function fig = slam_plot()

theta = linspace(pi,2*pi,200);
y_axis_span = linspace(0,1,120);
two_Tf_over_H_max = 3;
[THETA,Y_AXIS_SPAN] = meshgrid(theta,y_axis_span);

% this makes it so y_axis_span=0 corresponds to y=|sin(theta)| and
% y_axis_span=1 corresponds to y=two_Tf_over_H_max, to get a good sin shape
TWO_TF_OVER_H = abs(sin(THETA)) + Y_AXIS_SPAN .* (two_Tf_over_H_max - abs(sin(THETA)));

sqrt_term = sqrt(TWO_TF_OVER_H - sin(THETA).^2);
X_max_over_amp = sqrt_term + cos(THETA);
X_min_over_amp = -sqrt_term + cos(THETA);
idx_imag = sqrt_term~=real(sqrt_term);
X_max_over_amp(idx_imag) = NaN;
X_min_over_amp(idx_imag) = NaN;

fig = figure;

% x max plot
ax1 = subplot(1,2,1);
contours = sort([0 .1 -.1 linspace( min(X_max_over_amp(:)), max(X_max_over_amp(:)), 10 )]);
contourf(THETA/pi,TWO_TF_OVER_H,X_max_over_amp,contours)
colorbar;
xlabel('$\theta/\pi$','Interpreter','latex')
ylabel('$\displaystyle \frac{\Delta z_{slam}}{H/2}$','Interpreter','latex')
my_text = {'$\textrm{Infeasible if}~ \frac{\Delta z_{slam}}{H/2} < 1$','$\textrm{       and } |\theta-\pi|<\pi/2$'};
text(1.28,.4,my_text,'Interpreter','latex','FontSize',16)
colormap(ax1,bluewhitered)

p1 = [1.33,0.55];
p2 = [1.2,0.7];
arrow(p1,p2)

% x min plot
subplot 122
contours = sort([0 .1 -.1 linspace( min(X_min_over_amp(:)), max(X_min_over_amp(:)), 10 )]);
contourf(THETA/pi,TWO_TF_OVER_H,X_min_over_amp,contours)
colorbar;
xlabel('$\theta/\pi$','Interpreter','latex')
ylabel('$\displaystyle \frac{\Delta z_{slam}}{H/2}$','Interpreter','latex')
my_text = '$\textrm{Undefined if}~ \frac{\Delta z_{slam}}{H/2} < |\sin\theta|$';
text(1.19,.4,my_text,'Interpreter','latex','FontSize',16)

% plot aesthetics for both
improvePlot
set(gcf,"Position",[20 50 1500 800])
title(ax1,'$\frac{\xi_{max,slam}}{H/2}$','Interpreter','latex','FontSize',24)
title('$\frac{\xi_{min,slam}}{H/2}$','Interpreter','latex','FontSize',24)
colormap(gca,bluewhitered)

end
