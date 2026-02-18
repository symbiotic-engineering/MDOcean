function f = make_slamming_plot(T,Hs,phase_X,theta_slam,X_slam_max,X_slam_min,X_mag)

    f = figure;
    t = tiledlayout(1,5,'TileSpacing','tight');

    % map 0-2pi onto 0-pi without changing the value of cos(theta)
    % toy code to show this identity works:
    % x=linspace(0,2*pi,100); y=cos(x); x2=pi-abs(pi-x); y2=cos(x2); figure; plot(x,y,x,y2,'--'); figure; plot(x,x2)
    theta_new = pi-abs(pi-wrapTo2Pi(theta_slam));

    min_level = min([phase_X(:);theta_new(:)]);
    max_level = max([phase_X(:);theta_new(:)]);
    lims = [min_level,max_level]/pi;
    levels = sort([0.5 0.55 floor(lims(1)):0.1:ceil(lims(2))]);

    ax1 = nexttile;
    [C,h]=contourf(T,Hs,phase_X/pi,levels);
    clim(lims)
    xlabel('T_e'); 
    ylabel('H_s'); 
    title('$\angle\hat{\xi}/\pi$','Interpreter','latex'); 
    clabel(C,h);
    colormap(ax1,bluewhitered(256,.5))
    colorbar; 

    ax2 = nexttile;
    [C,h]=contourf(T,Hs,theta_new/pi,levels);
    clim(lims)
    xlabel('T_e'); 
    %ylabel('H_s'); 
    set(gca,'Yticklabel',[]) 
    title('$\theta/\pi$','Interpreter','latex'); 
    clabel(C,h);
    colormap(ax2,bluewhitered(256,.5))
    colorbar; 

    % min amplitude
    ax3 = nexttile;
    min_level = min(X_slam_min(:));
    max_level = max(X_slam_min(:));
    lims = [min_level,max_level];
    levels = sort([0 floor(lims(1)):0.5:-2.5 -2.2 -2.1 -2:0.2:ceil(lims(2))]);
    [C,h]=contourf(T,Hs,X_slam_min,levels);
    xlabel('T_e'); 
    %ylabel('H_s'); 
    set(gca,'Yticklabel',[]) 
    title('$\xi_{min,slam}$','Interpreter','latex');
    clabel(C,h);
    colormap(ax3,bluewhitered)
    colorbar
    hold on

    % max amplitude
    min_level = min([X_slam_max(:);X_mag(:)]);
    max_level = max([X_slam_max(:);X_mag(:)]);
    lims = [min_level,max_level];
    lims_nearest_half = fix(2*lims)/2;
    levels = sort([0 lims lims_nearest_half(1):0.5:lims_nearest_half(2)]);

    ax4 = nexttile;
    [C,h]=contourf(T,Hs,X_slam_max,levels);
    clim([min_level max_level])
    xlabel('T_e');
    %ylabel('H_s'); 
    set(gca,'Yticklabel',[]) 
    title('$\xi_{max,slam}$','Interpreter','latex'); 
    colorbar; 
    colormap(ax4,bluewhitered)
    clabel(C,h); 
    hold on

    % amplitude plot
    ax5 = nexttile;
    [C,h]=contourf(T,Hs,X_mag,levels);
    clim([min_level max_level])
    xlabel('T_e'); 
    %ylabel('H_s');
    set(gca,'Yticklabel',[]) 
    title('$|\hat{\xi}|$','Interpreter','latex'); 
    colorbar; 
    colormap(ax5,bluewhitered)
    clabel(C,h);
    hold on

    improvePlot;

    % hatches should happen after improveplot
    violating_min = X_slam_min>X_mag;
    nexttile(3)
    hatch(violating_min,ax3,T,Hs)
    
    violating_max = X_slam_max<X_mag;
    nexttile(4)
    hatch(violating_max,ax4,T,Hs)

    nexttile(5)
    hatch(violating_min|violating_max,ax5,T,Hs);

    f.Position(3) = 1500;
    title(t,'Slamming Model for Nominal Float Design')

end

function hatch(bool,ax,T,Hs)
     if any(bool(:))
        hold on
        tmp_clim = clim;
        [~,h2] = contourf(ax,T,Hs,bool,[1 1],'Fill','off');
        hh = hatchfill2(h2,'cross');
        hh.Color = [0 0 0 .5];
        clim(tmp_clim)
    end

end