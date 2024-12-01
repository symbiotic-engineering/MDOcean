function pareto_curve_heuristics()
    p = parameters();
    b = var_bounds();
    p_w = parameters('wecsim');
    b_w = var_bounds('wecsim');
    
    %[x,fval] = pareto_search();
    d=dir("**/pareto_search_results*");
    load(d(end).name)
    cols = b.idxs_recover;
    X = x(:,cols); % swap indices based on solver generated function
    X = [X ones(length(X),1)]; % add extra column for material 
    LCOE = fval(:,1);
    Pvar = fval(:,2);
    
    [LCOE, minLCOE, idx_best_LCOE, LCOE_nom,  LCOE_nom_sim,  LCOE_solar,  LCOE_balanced,...
     Pvar, minPvar, idx_best_Pvar, P_var_nom, P_var_nom_sim, P_var_solar, P_var_balanced,...
     x_best_LCOE, x_best_Pvar, x_nom, x_balanced, idxo] = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w);
    
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
                           Pvar, minPvar, idx_best_Pvar, X, idxo, p.LCOE_max, b.var_names_pretty(1:end-1))

end
%%
function [LCOE, minLCOE, idx_best_LCOE, LCOE_nom, ...
         LCOE_nom_sim,  LCOE_solar,  LCOE_balanced,...
         Pvar, minPvar, idx_best_Pvar, P_var_nom, ...
         P_var_nom_sim, P_var_solar, P_var_balanced,...
         x_best_LCOE, x_best_Pvar, x_nom, x_balanced, idxo] ...
                                        = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w)

    % find basic pareto front
    [~,idxo] = paretoFront(-1*[LCOE Pvar]);
    [minLCOE,idx_best_LCOE] = min(LCOE);
    [minPvar,idx_best_Pvar] = min(Pvar);
    
    % RM3 nominal - actual
    RM3_report = validation_inputs('report');
    RM3_wecsim = validation_inputs('wecsim');
    LCOE_nom = [RM3_report.LCOE(4) RM3_wecsim.LCOE(4)];
    P_var_nom = [RM3_report.c_v    RM3_wecsim.c_v];
    
    % RM3 nominal - simulated
    x_nom = [b.X_noms b_w.X_noms];
    x_nom(end+1,:) = 1;

    [LCOE_report_sim,P_var_report_sim] = simulation(x_nom(:,1),p);
    [LCOE_wecsim_sim,P_var_wecsim_sim] = simulation(x_nom(:,2),p_w);
    
    LCOE_nom_sim =  [LCOE_report_sim, LCOE_wecsim_sim];
    P_var_nom_sim = [P_var_report_sim,P_var_wecsim_sim];

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
    
    % RM3 nominal reference - report
    plot(LCOE_nom(1),P_var_nom(1),'rd','HandleVisibility','off')
    plot(LCOE_nom_sim(1),P_var_nom_sim(1),'rs','HandleVisibility','off')

    % RM3 nominal reference - wecsim
    plot(LCOE_nom(2),P_var_nom(2),'md','HandleVisibility','off')
    plot(LCOE_nom_sim(2),P_var_nom_sim(2),'ms','HandleVisibility','off')
    
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
    yellow = [1, .87, .2];
    plot(LCOE_solar, P_var_solar,'o','MarkerSize',12,'MarkerEdgeColor',yellow,...
        'MarkerFaceColor',yellow,'HandleVisibility','off')
    % for the yellow color to work, do not use improvePlot below here
    
    % text labels
    sz = 14;
    text(LCOE_nom(1)+.03,P_var_nom(1),'RM3 Report','FontSize',sz)
    text(LCOE_nom(1)+.01,P_var_nom(1)-5,'Actual [10]','FontSize',sz)
    text(LCOE_nom_sim(1)-.02,P_var_nom_sim(1)+5,'RM3 Report','FontSize',sz)
    text(LCOE_nom_sim(1)+.03,P_var_nom_sim(1),'Sim','FontSize',sz)

    text(LCOE_nom(2)+.03,P_var_nom(2),'RM3 WecSim','FontSize',sz)
    text(LCOE_nom(2)+.01,P_var_nom(2)-5,'Actual','FontSize',sz)
    text(LCOE_nom_sim(2)-.02,P_var_nom_sim(2)+5,'RM3 WecSim','FontSize',sz)
    text(LCOE_nom_sim(2)+.03,P_var_nom_sim(2),'Sim','FontSize',sz)

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
        
        % RM3 report
        axes('Position',[.7 .53 mini_plot_size])
        box on
        visualize_geometry(x_nom(:,1),p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])

        % RM3 WecSim
        axes('Position',[.8 .7 mini_plot_size])
        box on
        visualize_geometry(x_nom(:,2),p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
    end
end

%%
function [] = design_heuristics_plot(overallLCOE, minLCOE, idx_best_LCOE, x_best_LCOE, ...
                                     overallPvar, minPvar, idx_best_Pvar, ...
                                     overallX, idxo, LCOE_max, var_names)

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
    
    X_pareto_sorted_scaled = X_pareto_sorted_scaled(:,1:end-1); % get rid of material
    
    windowSize = round(length(idx_sort) * 5/100);
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    y = zeros(size(X_pareto_sorted_scaled));
    for i=1:size(X_pareto_sorted_scaled,2)
        x = X_pareto_sorted_scaled(:,i); 
        x_padded = [ones(windowSize,1); x];
        yy = filter(b,a,x_padded); % moving average filter
        y(:,i) = yy((windowSize+1):end);
    end
    
    % unfiltered
    figure
    semilogy(pct_angle,X_pareto_sorted_scaled)
    title('Unfiltered Design Heuristics')
    xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Optimal Design Value')
    legend(var_names,'Location','eastoutside')
    ylim([0 15])
    improvePlot
    grid on
    set(gca,'YMinorGrid','on')
    
    % filtered
    cols = {'r:','r--','r-','r-.','b:','b--','b-','g:','g--','g-','g-.','g.','g*','b.'};
    figure
    for i=1:size(X_pareto_sorted_scaled,2)
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
    legend(var_names,'Location','eastoutside')
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
    improvePlot
    cent = char(0162);
    legend(['LCOE (' cent '/kWh)'],'c_v (%)',Location='northeast')
end