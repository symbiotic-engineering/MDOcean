function make_slamming_plot(T,Hs,theta_slam,X_below_wave,X_star,X_below_simple)

    f = figure;
    t = tiledlayout(1,4);

    nexttile
    contourf(T,Hs,theta_slam/pi);
    xlabel('T_e'); 
    ylabel('H_s'); 
    title('\theta/\pi'); 
    colorbar; 
    improvePlot;

    ax = nexttile;
    contourf(T,Hs,X_star,20);
    xlabel('T_e'); 
    ylabel('H_s'); 
    title('X^*'); 
    improvePlot; 
    colormap(ax,bluewhitered)
    colorbar

    levels = [-.1 0 .05 .25 .5 1:5 10 20 50 100];

    ax = nexttile;
    [C,h]=contourf(T,Hs,X_below_wave,levels);
    clim([-.1 5]);
    xlabel('T_e');
    ylabel('H_s'); 
    title('X_{slam}/X - 1'); 
    colorbar; 
    colormap(ax,bluewhitered)
    clabel(C,h); 
    improvePlot;

    ax = nexttile;
    [C,h]=contourf(T,Hs,X_below_simple,levels);
    clim([-.1 5]); 
    xlabel('T_e'); 
    ylabel('H_s');
    title('(T_f-H/2)/X - 1'); 
    colorbar; 
    colormap(ax,bluewhitered)
    clabel(C,h);
    improvePlot;     

    f.Position(3) = 1500;
    title(t,'Slamming Model for Min LCOE Design')
end
