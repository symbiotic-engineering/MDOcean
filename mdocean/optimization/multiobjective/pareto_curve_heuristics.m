function pareto_curve_heuristics()
    p = parameters();
    b = var_bounds();
    
    %[x,fval] = pareto_search();
    load("pareto_search_results.mat")
    cols = [1 3 6 5 4 2 7];
    X = x(:,cols); % swap indices based on solver generated function
    X = [X ones(length(X),1)]; % add eigth column for material 
    LCOE = fval(:,1);
    Pvar = fval(:,2);
    
    [LCOE, minLCOE, idx_best_LCOE, LCOE_nom,  LCOE_nom_sim,  LCOE_solar,  LCOE_balanced,...
     Pvar, minPvar, idx_best_Pvar, P_var_nom, P_var_nom_sim, P_var_solar, P_var_balanced,...
     x_best_LCOE, x_best_Pvar, x_nom, x_balanced, idxo] = process_pareto_front(LCOE,Pvar,X,p,b);
    
    %% super simple "pareto" plot of just single objective optimizations
    showSingleObj = true;
    showImages = false;
    pareto_plot(LCOE, minLCOE, idx_best_LCOE, LCOE_nom,  LCOE_nom_sim,  LCOE_solar,  NaN*LCOE_balanced,...
                Pvar, minPvar, idx_best_Pvar, P_var_nom, P_var_nom_sim, P_var_solar, NaN*P_var_balanced,...
                x_best_LCOE, x_best_Pvar, x_nom, x_balanced, [], showSingleObj, showImages, p)
    
    %% simple pareto plot
    showSingleObj = false;
    showImages = false;
    pareto_plot(LCOE, minLCOE, idx_best_LCOE, LCOE_nom,  LCOE_nom_sim,  LCOE_solar,  LCOE_balanced,...
                Pvar, minPvar, idx_best_Pvar, P_var_nom, P_var_nom_sim, P_var_solar, P_var_balanced,...
                x_best_LCOE, x_best_Pvar, x_nom, x_balanced, idxo, showSingleObj, showImages, p)
    
    %% plot pareto front with annotations and embedded images of three recommended designs
    showSingleObj = true;
    showImages = true;
    pareto_plot(LCOE, minLCOE, idx_best_LCOE, LCOE_nom,  LCOE_nom_sim,  LCOE_solar,  LCOE_balanced,...
                Pvar, minPvar, idx_best_Pvar, P_var_nom, P_var_nom_sim, P_var_solar, P_var_balanced,...
                x_best_LCOE, x_best_Pvar, x_nom, x_balanced, idxo, showSingleObj, showImages, p)
    
    %% plots for DVs as a fn of percent along the pareto
    design_heuristics_plot(LCOE, minLCOE, idx_best_LCOE, x_best_LCOE, ...
                           Pvar, minPvar, idx_best_Pvar, X, idxo, p.LCOE_max)

end
%%
function [LCOE, minLCOE, idx_best_LCOE, LCOE_nom, ...
         LCOE_nom_sim,  LCOE_solar,  LCOE_balanced,...
         Pvar, minPvar, idx_best_Pvar, P_var_nom, ...
         P_var_nom_sim, P_var_solar, P_var_balanced,...
         x_best_LCOE, x_best_Pvar, x_nom, x_balanced, idxo] ...
                                        = process_pareto_front(LCOE,Pvar,X,p,b)

    % find basic pareto front
    [~,idxo] = paretoFront(-1*[LCOE Pvar]);
    [minLCOE,idx_best_LCOE] = min(LCOE);
    [minPvar,idx_best_Pvar] = min(Pvar);
    
    % RM3 nominal - actual
    RM3 = validation_inputs();
    LCOE_nom = RM3.LCOE(4);
    P_var_nom = RM3.c_v;
    
    % RM3 nominal - simulated
    x_nom = b.X_noms;
    x_nom(8) = 1;
    p_nom = p;
    %p_nom.power_max = RM3.power_max;
    [LCOE_nom_sim,P_var_nom_sim] = simulation(x_nom,p_nom);
    
    % solar
    LCOE_solar = 0.03;
    P_var_solar = 125;
    
    % balanced design
    [~,idx_balanced] = min(abs(Pvar-100));
    LCOE_balanced = LCOE(idx_balanced);
    P_var_balanced = Pvar(idx_balanced);
    
    % idenitfy design variables for best designs
    x_best_LCOE = X(idx_best_LCOE,:)
    x_best_Pvar = X(idx_best_Pvar,:)
    x_balanced = X(idx_balanced,:)
end

%%
function [] = pareto_plot(LCOE,minLCOE,idx_best_LCOE,LCOE_nom, LCOE_nom_sim, LCOE_solar, LCOE_balanced,...
                          Pvar,minPvar,idx_best_Pvar,P_var_nom,P_var_nom_sim,P_var_solar,P_var_balanced,...
                          x_best_LCOE,x_best_Pvar,x_nom,x_balanced,idxo,showSingleObj,showImages,p)
    figure
    % overall pareto front
    plot(LCOE(idxo),Pvar(idxo),'bs','MarkerFaceColor','b','HandleVisibility','off')
    hold on
    
    % utopia point
    plot(minLCOE,minPvar,'gp','MarkerFaceColor','g','MarkerSize',20,'HandleVisibility','off')
    
    % RM3 nominal reference
    plot(LCOE_nom,P_var_nom,'rd','HandleVisibility','off')
    plot(LCOE_nom_sim,P_var_nom_sim,'rs','HandleVisibility','off')
    
    if showSingleObj  
        % black squares for 3 ref points
        plot(minLCOE,Pvar(idx_best_LCOE),'ks','HandleVisibility','off')
        plot(LCOE(idx_best_Pvar),minPvar,'ks','HandleVisibility','off')
        plot(LCOE_balanced,P_var_balanced,'ks','HandleVisibility','off')
    end
    
    % axis labels
    xlabel('LCOE ($/kWh)')
    ylabel('Power Variation (%)')
    title('Pareto Front')
    xlim([0 1])
    ylim([25 165])
    improvePlot
    
    % solar reference
    plot(LCOE_solar, P_var_solar,'o','MarkerSize',12,'MarkerEdgeColor',[1, .87, .2],'MarkerFaceColor',[1, .87, .2],'HandleVisibility','off')
    % for the yellow color to work, do not use improvePlot below here
    
    % text labels
    sz = 14;
    text(LCOE_nom+.03,P_var_nom,'Nominal','FontSize',sz)
    text(LCOE_nom+.01,P_var_nom-5,'Actual [10]','FontSize',sz)
    text(LCOE_nom_sim-.02,P_var_nom_sim+5,'Nominal','FontSize',sz)
    text(LCOE_nom_sim+.03,P_var_nom_sim,'Sim','FontSize',sz)
    text(minLCOE+.03,minPvar-2,'Utopia Point','FontSize',sz)
    text(LCOE_solar+.03,P_var_solar,'Solar','FontSize',sz)
    if showSingleObj
        text(LCOE(idx_best_LCOE)+.03,Pvar(idx_best_LCOE),'Cheapest','FontSize',sz)
        text(LCOE(idx_best_Pvar)+.03,Pvar(idx_best_Pvar)-3,'Least Variable','FontSize',sz)
        text(LCOE_balanced-.15,P_var_balanced+5,'Balanced Design','FontSize',sz)
    end

    if showImages
        mini_plot_size = [.2 .22];
        % small corner pictures of best geometries
        % upper left
        axes('Position',[.28 .6 mini_plot_size])
        box on
        visualize_geometry(x_best_LCOE,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
        
        % lower right
        axes('Position',[.51 .23 mini_plot_size])
        box on
        visualize_geometry(x_best_Pvar,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
        
        % balanced
        axes('Position',[.10 .28 mini_plot_size])
        box on
        visualize_geometry(x_balanced,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
        
        % RM3
        axes('Position',[.7 .53 mini_plot_size])
        box on
        visualize_geometry(x_nom,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
    end
end

%%
function [] = design_heuristics_plot(overallLCOE, minLCOE, idx_best_LCOE, x_best_LCOE, ...
                                     overallPvar, minPvar, idx_best_Pvar, ...
                                     overallX, idxo, LCOE_max)

    LCOE_pareto = overallLCOE(idxo);
    Pvar_pareto = overallPvar(idxo);
    [LCOE_pareto_sorted,idx_sort] = sort(LCOE_pareto(LCOE_pareto<LCOE_max));
    Pvar_pareto_sorted = Pvar_pareto(LCOE_pareto<LCOE_max);
    Pvar_pareto_sorted = Pvar_pareto_sorted(idx_sort);
    
    pct = linspace(0,100,length(idx_sort));
    %pct_angle = 100/(pi/2) * atan((LCOE_pareto_sorted - minLCOE) ./ (Pvar_pareto_sorted - minPvar));
    num = (overallPvar(idx_best_LCOE) - Pvar_pareto_sorted) / (overallPvar(idx_best_LCOE) - minPvar);
    den = (overallLCOE(idx_best_Pvar) - LCOE_pareto_sorted) / (overallLCOE(idx_best_Pvar) - minLCOE);
    pct_angle = 100/(pi/2) * atan( num ./ den);
    pct_angle(pct_angle==-100) = 100;
    
    X_pareto = overallX(idxo,:);
    X_pareto = X_pareto(LCOE_pareto<LCOE_max,:);
    X_pareto_sorted = X_pareto(idx_sort,:);
    
    X_pareto_sorted_scaled = X_pareto_sorted ./ repmat(x_best_LCOE,length(idx_sort),1);
    
    X_pareto_sorted_scaled = X_pareto_sorted_scaled(:,1:7); % get rid of material
    
    windowSize = round(length(idx_sort) * 5/100);
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    y = zeros(size(X_pareto_sorted_scaled));
    for i=1:7
        x = X_pareto_sorted_scaled(:,i); 
        x_padded = [ones(windowSize,1); x];
        yy = filter(b,a,x_padded); % moving average filter
        y(:,i) = yy((windowSize+1):end);
    end
    
    var_names_pretty = {'D_f',...    % outer diameter of float (m)	
                'D_s/D_f',...
                'h_f/D_f',...      	
                'T_s/h_s',...      	
                'F_{max}',...     % max force (N)
                'B_p',...     % internal damping of controller (Ns/m)	
                'w_n'};         % natural frequency (rad/s)
    
    
    % unfiltered
    figure
    semilogy(pct_angle,X_pareto_sorted_scaled)
    title('Unfiltered Design Heuristics')
    xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Optimal Design Value')
    legend(var_names_pretty,'Location','eastoutside')
    ylim([0 15])
    improvePlot
    grid on
    set(gca,'YMinorGrid','on')
    
    % filtered
    cols = {'r:','r--','r-','r-.','b:','b--','b-'};
    figure
    for i=1:7
        semilogy(pct_angle,y(:,i),cols{i})
        hold on
    end
    
    % make fake major grid lines (didn't do grid on because then major lines show up for .2 .5 2 5 too)
    x_grid = [3 97];
    y_grid = [.1 1];
    y_tick = [.01 .02 .05 .1 .2 .5 1 2 5];
    for i=1:length(y_grid)
        plot(x_grid, y_grid(i)*[1 1],'Color',[.85 .85 .85]);
    end
    
    title('Design Heuristics')
    %xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Optimal Design Value')
    ylim([.01 5])
    set(gca,'YTick',y_tick)
    improvePlot
    legend(var_names_pretty,'Location','eastoutside')
    set(gca,'YMinorGrid','on')
    set(gca,'XGrid','on')
    set(gca, 'Children', flipud(get(gca, 'Children')) ) % put fake gridlines behind real lines
    
    figure
    plot(pct_angle,LCOE_pareto_sorted*100,'Color',[0.4940 0.1840 0.5560]) % purple
    hold on
    plot(pct_angle,Pvar_pareto_sorted,'Color',[0.4660 0.6740 0.1880]) % green
    grid on
    xlabel('Percent along the Pareto Curve')
    ylabel('Objective Value')
    pct = [0 4.6683 16.1186 27.2305 28.6278 30.4717 34.0738 38.4791 40.3735 41.0267 ...
    41.9764 43.3829 43.5881 44.4236 44.8667 45.7298 50.1799 51.8743 52.9700 55.1960 ...
    56.1828 58.0417 58.8600 60.9711 61.2873 61.5743 62.1110 62.7047 63.3152 64.3656 ...
    65.9389 66.9220 67.6261 68.5065 68.7219 69.0767 70.5319 70.8866 71.3303 72.1607 ...
    72.2847 73.9855 75.2496 78.0466 78.7258 80.3724 83.5218 84.3179 84.6826 86.4844 ...
    87.0652 88.0020 89.1332 89.8590 91.8756 96.5597 97.7855 99.0541]';
    lcoe = [10.0000 14.2494 14.9043 16.3479 16.5465 16.8356 17.5015 18.5171 19.0298 ...
    19.2161 19.4955 19.9336 20.0000 20.2759 20.4259 20.7251 22.4162 23.1281 23.6075 ...
    24.6018 25.0557 25.9355 26.3336 27.3902 27.5518 27.6990 27.9760 28.2847 28.6044 ...
    29.1591 30.0000 30.5316 30.9152 31.3981 31.5409 31.7360 32.5223 32.7232 32.9696 ...
    33.4381 33.5070 34.4715 35.1936 36.8083 37.2151 38.1557 40.0000 40.4700 40.6851 ...
    41.7548 42.1005 42.6605 43.3483 43.7747 44.9983 47.8588 48.6169 49.4073]';
    pvar = [158.8509 151.7591 134.3300 117.4241 115.2494 112.3817 106.8146 100.1354 ...
    97.3378 96.3835 95.0058 93.0008 92.7122 91.5459 90.9332 89.7518 83.9438 81.8848 ...
    80.6024 78.0478 76.9518 74.9637 74.1240 72.0641 71.7683 71.5025 71.0121 70.4791 ...
    69.9409 69.0366 67.7310 66.9463 66.3994 65.7339 65.6954 65.4326 64.2794 64.0486 ...
    63.7328 63.1846 63.0969 62.0126 61.2449 59.6798 59.4164 58.4211 56.8669 56.5049 ...
    56.3351 55.5597 55.3120 54.9383 54.6481 54.2157 53.5928 51.9609 51.5273 51.1327]';
    plot(pct,lcoe,'Color',[1 0 0])
    plot(pct,pvar,'Color',[0 0 1])
    improvePlot
    cent = char(0162);
    legend(['LCOE (' cent '/kWh)'],'c_v (%)', ...
        'LCOE, const geom','Pvar, const geom',Location='northeast')
end