function fig = slam_plot()

theta = linspace(0,pi,200);
T_f_over_H_span = linspace(0,1,120);
two_Tf_over_H_max = 3;
[THETA,TWO_TF_OVER_H_SPAN] = meshgrid(theta,T_f_over_H_span);

TWO_TF_OVER_H = sin(THETA) + TWO_TF_OVER_H_SPAN .* (two_Tf_over_H_max - sin(THETA));
squared = sin(THETA)./(TWO_TF_OVER_H);
X_map = TWO_TF_OVER_H .* (sqrt(1 - squared.^2) - 1) - cos(THETA);
idx_imag = X_map~=real(X_map);
X_map(idx_imag) = NaN;%-2*TF_OVER_H(idx_imag);

fig = figure;
contourf(THETA/pi,TWO_TF_OVER_H,real(X_map),10)
cb=colorbar;
cb.Label.String = '$X^*$';
cb.Label.Interpreter = 'latex';
cb.Label.Rotation = 0;
cb.Label.FontSize = 24;
xlabel('$\theta/\pi$','Interpreter','latex')
ylabel('$\frac{T_f}{H/2}$','Interpreter','latex')
clim([-1 1])
%title()
my_text = '$X_{slam}=0 ~\textrm{if}~ \frac{T_f}{H/2} < \sin\theta$';
text(.19,.4,my_text,'Interpreter','latex','FontSize',24)
improvePlot
set(gcf,"Position",[100 100 750 600])
title('$X_{slam}=T_f+\frac{H}{2} X^*$','Interpreter','latex','FontSize',24)

end