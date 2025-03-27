function plot_power_matrix(X,p,filename_uuid)

[CW_over_CW_max, P_wave, CW_max, P_elec] = check_max_CW(filename_uuid, X, p);

[T,Hs] = meshgrid(p.T,p.Hs);
figure

subplot(2,5,1)
sub_one = subplot(251);
sub_one.Position = [0.05, 0.5838, 0.24, 0.3412];
contourf(T,Hs,P_wave);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Raw Wave Power Density (W/m)', 'FontSize', 12, 'Position', [11, 7, 0])
colorbar

subplot(2,5,2)
sub_two = subplot(252);
sub_two.Position = [0.315, 0.5838, 0.01, 0.3412];
text(0.5,0.5,'x','FontSize',30)
axis off

subplot(2,5,3)
sub_three = subplot(253);
sub_three.Position = [0.3900, 0.5838, 0.24, 0.3412];
contourf(T,Hs,CW_over_CW_max)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Device Capture Efficiency','FontSize',12, 'Position', [11.75, 7, 0])
colorbar

subplot(2,5,4)
sub_four = subplot(254);
sub_four.Position = [0.650, 0.5838, 0.01, 0.3412];
text(0.5,0.5,'x','FontSize',30)
axis off

subplot(2,5,5)
sub_five = subplot(255);
sub_five.Position = [0.7250, 0.5838, 0.23, 0.3412];
contourf(T,Hs,CW_max)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Radiation Capture Width Limit (m)','FontSize',12, 'Position', [11.75, 7, 0])
colorbar

subplot(2,5,6)
sub_six = subplot(256);
sub_six.Position = [0.19, 0.1100, 0.01, 0.3412];
text(0.55,0.5,'x','FontSize',30)
axis off

subplot(2,5,7)
sub_seven = subplot(257);
sub_seven.Position = [0.2628, 0.1100, 0.26, 0.3412];
contourf(T,Hs,p.JPD)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Probability at Site (%)','FontSize',12)
colorbar

subplot(2,5,8)
sub_eight = subplot(258);
sub_eight.Position = [0.5496, 0.1100, 0.01, 0.3412];
text(0.52,0.5,'=','FontSize',30)
axis off

subplot(2,5,9)
sub_nine=subplot(259);
sub_nine.Position=[0.6350, 0.1100, 0.26, 0.3412];
P_product = P_wave .* CW_over_CW_max .* CW_max .* p.JPD;
contourf(sub_nine,T,Hs,P_product)
eff = P_elec ./ P_product;
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Weighted Power (W)','FontSize',12)
colorbar

improvePlot

end
