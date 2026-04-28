function [carbon_contour,price_contour] = timeseries_to_sea_state_matrix(Hs,T,jpd_Hs,jpd_T,carbon_interpolated,price_interpolated,ploton)

%make the boundary for the first bucket 0
p.Hs = [0; jpd_Hs];
p.T = [0, jpd_T];


%find the average price or carbon value for each bin
price_contour=zeros(length(p.Hs)-1,length(p.T)-1);
carbon_contour=zeros(length(p.Hs)-1,length(p.T)-1);
count = 1;
for i = 1:length(p.Hs)-1
    for j = 1:length(p.T)-1
        idx_use = Hs > p.Hs(i) & Hs < p.Hs(i+1) & T > p.T(j) & T < p.T(j+1);
        price_contour(i,j) = mean(price_interpolated(idx_use),'all');
        carbon_contour(i,j) = mean(carbon_interpolated(idx_use),'all');
        price_var(count) = var(price_interpolated(idx_use));
        carbon_var(count) = var(carbon_interpolated(idx_use));
        count = count + 1;
    end
end

if ploton == true
    %graph the bins
    figure
    hold on
    for i = 1:length(p.Hs)
         yline(p.Hs(i))
    end
    for j = 1:length(p.T)
         xline(p.T(j))
    end
    xlabel("Wave Period (s)")
    ylabel("Wave Height (m)")
        title('Wave Height vs. Wave Period')
    xticks(p.T)
    yticks(p.Hs)
    hold off
    improvePlot

    %plot timeseries
    t = linspace(1,8760,8760);
    subplot(2,2,1)
    plot(t, T)
    xlabel('Time (hr)')
    ylabel('Wave Period (s)')
    title('Wave Period Time Series')
    xlim([0,8760])
    improvePlot

    subplot(2,2,2)
    plot(t, Hs)
    xlabel('Time (hr)')
    ylabel('Wave Height (m)')
    title('Wave Height Time Series')
    xlim([0,8760])
    improvePlot

    subplot(2,2,3)
    plot(t, price_interpolated)
    xlabel('Time (hr)')
    ylabel('Energy Price ($/kWh)')
    title('Energy Price Time Series')
    xlim([0,8760])
    improvePlot

    subplot(2,2,4)
    plot(t, carbon_interpolated)
    xlabel('Time (hr)')
    ylabel('Emitted Carbon (kg/kWh)')
    title('Emitted Carbon Time Series')
    xlim([0,8760])
    improvePlot

    %make contour plots
    [T_mdocean,Hs_mdocean] = meshgrid(p.T(2:end),p.Hs(2:end));
    figure
    subplot(1,2,1)
    contourf(T_mdocean, Hs_mdocean, price_contour)
    xlabel('T (s)')
    ylabel('H_{s} (m)')
    title('Price ($/kWh)')
    colorbar
    improvePlot
    
    subplot(1,2,2)
    contourf(T_mdocean, Hs_mdocean, carbon_contour)
    title('Carbon (kg/kWh)')
    xlabel('T (s)')
    ylabel('H_{s} (m)')
    colorbar
    improvePlot
end