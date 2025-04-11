function plot_power_matrix(X,p,filename_uuid)

[CW_over_CW_max, P_wave, CW_max, P_elec] = check_max_CW(filename_uuid, X, p);

[T,Hs] = meshgrid(p.T,p.Hs);
figure

subplot(2,6,1)
sub_one = subplot(261);
sub_one.Position = [0.06, 0.5838, 0.12, 0.3412];
contourf(sub_one,T,Hs,P_wave);
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Raw Wave Power Density (W/m)', 'FontSize', 10, 'Position', [11, 7, 0])
colorbar

subplot(2,6,2)
sub_two = subplot(262);
sub_two.Position = [0.24, 0.5838, 0.01, 0.3412];
text(0.5,0.5,'x','FontSize',40)
axis off

subplot(2,6,3)
sub_three = subplot(263);
sub_three.Position = [0.36, 0.5838, 0.12, 0.3412];
contourf(sub_three,T,Hs,CW_over_CW_max)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Device Capture Efficiency','FontSize',10, 'Position', [11.75, 7, 0])
colorbar

subplot(2,6,4)
sub_four = subplot(264);
sub_four.Position = [0.54, 0.5838, 0.01, 0.3412];
text(0.5,0.5,'x','FontSize',40)
axis off

subplot(2,6,5)
sub_five = subplot(265);
sub_five.Position = [0.68, 0.5838, 0.12, 0.3412];
contourf(sub_five,T,Hs,CW_max)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Radiation Capture Width Limit (m)','FontSize',10, 'Position', [11.75, 7, 0])
colorbar

subplot(2,6,6)
sub_six = subplot(266);
sub_six.Position = [0.92, 0.5838, 0.01, 0.3412];
text(0.55,0.5,'x','FontSize',40)
axis off

subplot(2,6,7)
sub_seven = subplot(267);
sub_seven.Position = [0.1, 0.1100, 0.12, 0.3412];
contourf(sub_seven,T,Hs,p.JPD)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Probability at Site (%)','FontSize',10)
colorbar

subplot(2,6,8)
sub_eight = subplot(268);
sub_eight.Position = [0.28, 0.1100, 0.01, 0.3412];
text(0.52,0.5, 'x','FontSize',40)
axis off

subplot(2,6,9)
sub_nine=subplot(269);
sub_nine.Position=[0.4269, 0.1100, 0.12, 0.3412];
P_product = P_wave .* CW_over_CW_max .* CW_max .* p.JPD;
eff = P_elec ./ P_product;
contourf(sub_nine,T,Hs,eff)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
title('Efficiency','FontSize',10)
colorbar

sub_ten=subplot(2,6,10);
sub_ten.Position = [0.6315, 0.1100, 0.01, 0.3412];
text(0.52,0.5, '=','FontSize',40)
axis off

sub_eleven=subplot(2,6,11);
sub_eleven.Position=[0.756, 0.1100, 0.12, 0.3412];
P_product_new = P_product .* eff;
contourf(sub_eleven,T,Hs,P_product_new)
xlabel('Wave Period T (s)')
ylabel('Wave Height Hs (m)')
colorbar
improvePlot

end
