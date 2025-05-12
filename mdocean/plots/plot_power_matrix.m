function plot_power_matrix(X,p,b,filename_uuid)

[CW_over_CW_max, P_wave, CW_max, P_elec, ...
    force_sat_ratio, drag_ratio, eff] = check_max_CW(filename_uuid, X, p, b, false);

[T,Hs] = meshgrid(p.T,p.Hs);
fig = figure;

P_elec_calc = P_wave .* CW_max .* CW_over_CW_max .* drag_ratio .* force_sat_ratio .* eff;
err = (P_elec - P_elec_calc) ./ P_elec;
assert(all(abs(err(~isnan(err))) < 1e-3,'all'))

t = tiledlayout(2,12);
t.TileSpacing = 'tight';
t.Padding = 'compact';
nexttile([1 2])
mycontour(T,Hs,P_wave/1000);
title('Wave Power (kW/m)', 'FontSize', 20, 'Position', [11, 7, 0])
c = colorbar;
% c.Position(1) = c.Position(1)-0.025;
% c.Position(3) = c.Position(3)/2.7;
% c.Position(4) = c.Position(4)/1.09;

nexttile;
text(0,0.5,'x','FontSize',40)
axis off

nexttile([1 2])
mycontour(T,Hs,CW_max)
title('Max Capture Width (m)','FontSize',20, 'Position', [11.75, 7, 0])
colorbar

nexttile
text(0,0.5,'x','FontSize',40)
axis off

nexttile([1 2])
mycontour(T,Hs,CW_over_CW_max*100)
title('Radiation Efficiency (%)','FontSize',20, 'Position', [11.75, 7, 0])
colorbar
caxis([0 100])

nexttile
text(0,0.5,'x','FontSize',40)
axis off

nexttile([1 2])
mycontour(T,Hs,drag_ratio*100)
title('Drag Efficiency (%)','FontSize',20, 'Position', [11.75, 7, 0])
colorbar
caxis([0 100])

nexttile
text(0,0.5,'x','FontSize',40)
axis off

nexttile([1 2])
mycontour(T,Hs,force_sat_ratio*100)
title('F_{max} Factor (%)','FontSize',20)
colorbar
caxis([0 100])

nexttile
text(0,0.5,'x','FontSize',40)
axis off

nexttile([1 2])
mycontour(T,Hs,eff*100)
title('Electrical Efficiency (%)','FontSize',20)
colorbar
caxis([0 100])

nexttile
text(0,0.5, 'x','FontSize',40)
axis off

nexttile([1 2])
p.JPD(p.JPD>0 & p.JPD < .001) = 0;
hrs_in_yr = 8766;
hours = p.JPD/100 * hrs_in_yr;
mycontour(T,Hs,p.JPD)
title('Site Probability (%)','FontSize',20)
colorbar

nexttile
text(0,0.5, '=','FontSize',40)
axis off

nexttile([1 2])
energy = P_elec .* hours;
mycontour(T,Hs,energy/1e6)
title('Annual Energy (MWh)')
colorbar

xlabel(t,'Wave Period T_e (s)')
ylabel(t,'Wave Height H_s (m)')
title(t,' ')
improvePlot

set(fig,'Position',[0 123 1530 620])

end


function mycontour(X,Y,Z)
    if numel(unique(Z(~isnan(Z)))) == 1
        % avoid "contour not rendered for constant zdata"
        x = [min(X,[],'all') max(X,[],'all')];
        y = [min(Y,[],'all') max(Y,[],'all')];
        imagesc('XData',x,'YData',y,'CData',Z,'AlphaData',~isnan(Z))
    else
        contourf(X,Y,Z)
    end

end

