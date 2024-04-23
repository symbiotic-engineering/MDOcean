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
    lambdaActive = lambda.active';
    lambdaLower = lambda.lower';
    lambdaUpper = lambda.upper';
    
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

    %% plots Lagrange multipliers as a fn of percent along the pareto
    lagrange_multiplier_plot(lambdaActive,lambdaUpper,lambdaLower)
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
    [~,idx_balanced] = min(abs(Pvar-40));
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
    plot(LCOE(idxo),Pvar(idxo),'bs','MarkerFaceColor','b')
    hold on
    
    % utopia point
    plot(minLCOE,minPvar,'gp','MarkerFaceColor','g','MarkerSize',20)
    
    % RM3 nominal reference
    plot(LCOE_nom,P_var_nom,'rd')
    plot(LCOE_nom_sim,P_var_nom_sim,'rs')
    
    if showSingleObj  
        % black squares for 3 ref points
        plot(minLCOE,Pvar(idx_best_LCOE),'ks')
        plot(LCOE(idx_best_Pvar),minPvar,'ks')
        plot(LCOE_balanced,P_var_balanced,'ks')
    end
    
    % axis labels
    xlabel('LCOE ($/kWh)')
    ylabel('Power Variation (%)')
    title('Pareto Front')
    xlim([0 1])
    ylim([25 165])
    improvePlot
    
    % solar reference
    plot(LCOE_solar, P_var_solar,'o','MarkerSize',12,'MarkerEdgeColor',[1, .87, .2],'MarkerFaceColor',[1, .87, .2])
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
        text(LCOE_balanced+.02,P_var_balanced+5,'Balanced Design','FontSize',sz)
    end
    
    if showImages
        mini_plot_size = [.2 .22];
        % small corner pictures of best geometries
        % upper left
        axes('Position',[.28 .7 mini_plot_size])
        box on
        visualize_geometry(x_best_LCOE,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
        
        % lower right
        axes('Position',[.52 .15 mini_plot_size])
        box on
        visualize_geometry(x_best_Pvar,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
        
        % balanced
        axes('Position',[.22 .25 mini_plot_size])
        box on
        visualize_geometry(x_balanced,p,true);
        set(gca,'XTickLabel',[],'YTickLabel',[])
        
        % RM3
        axes('Position',[.7 .5 mini_plot_size])
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
    improvePlot
    cent = char(0162);
    legend(['LCOE (' cent '/kWh)'],'c_v (%)',Location='northeast')
end

%%
function [] = lagrange_multiplier_plot(lambdaActive, ...
                                 lambdaUpper,lambdaLower)
    figure
    hold on
    [m,n]=size(lambdaActive);
    color_cell = {'.r','.g','.b','.c','.m','.y','.k',"#FF0000",'ob','+b','*b','xb','-b','|b'};
    %color_cell = {'r','g','b','c','m','y','k',[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    for i=1:m
        for j = 1:n
            %color = color_cell{rem(j,15)};
            plot(i*100/60,lambdaActive(i,j),'.','MarkerSize',12);
        end
    end
    %plot(pct_angle,lambda_active_sorted)
    title("Lagrange Multipliers",'FontSize',20)
    xlabel("Percent Along the Pareto Front",'FontSize',15)
    ylabel("Lagrange Multiplier",'FontSize',15)
    legend("prevent float too heavy","prevent float too light", ...
                     "prevent spar too heavy","prevent spar too light","stability",...
                     "float survives max force","spar survives max force",...
                     "damping plate survives max force","spar doesn't buckle", ...
                     'positive power','damping plate diameter','prevent float rising above top of spar',...
                     'prevent too expensive','prevent irrelevant max force')
    hold off

    % lower bound plot
    figure
    hold on
    [m,n]=size(lambdaLower);
    % I'll change these later
    for i=1:m
        for j = 1:n
            color = color_cell{rem(j,15)};
            plot(i*100/60,lambdaLower(i,j),color,'MarkerSize',12);
        end
    end
        title('Lower Bound Active Lagrange Multipliers')
        xlabel('Percent Along the Pareto Curve')
        ylabel('Lagrange Multiplier')
        xlim([-1,101]);
        ylim([-.05,.2])
        legend('WEC Surface Float Outer Diameter', ...
        'Ratio of WEC Surface Float Inner Diameter to Outer Diameter', ...
        'Ratio of WEC Surface Float Height to Outer Diameter', ...
        'Percent of WEC Spar Submergence','Maximum Powertrain Force', ...
        'Powertrain/Controller Damping','Controller Natural Frequency')
        hold off

    % upper bound plot
    figure
    hold on
    [m,n]=size(lambdaUpper)
    % I'll change these later
    for i=1:m
        for j = 1:n
            color = color_cell{rem(j,15)};
            plot(i*100/60,lambdaUpper(i,j),color,'MarkerSize',12);
        end
    end
        title('Upper Bound Active Lagrange Multipliers')
        xlabel('Percent Along the Pareto Curve')
        ylabel('Lagrange Multiplier')
        %xlim([-1,61]);
        ylim([-.02,.04])
        legend('WEC Surface Float Outer Diameter', ...
        'Ratio of WEC Surface Float Inner Diameter to Outer Diameter', ...
        'Ratio of WEC Surface Float Height to Outer Diameter', ...
        'Percent of WEC Spar Submergence','Maximum Powertrain Force', ...
        'Powertrain/Controller Damping','Controller Natural Frequency')
    improvePlot
    hold off
end