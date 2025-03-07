function pareto_curve_heuristics()
    p = parameters();
    b = var_bounds();
    p_w = parameters('wecsim');
    b_w = var_bounds('wecsim');
    
    %[x,fval] = pareto_search();
    d=dir("**/pareto_search_results*");
    load(d(end).name)

    if ~exist('tol','var')
        tol = 1e-6;
    end
    constraint_active_plot(residuals,fval,tol,b)

    cols = b.idxs_recover;
    X = x(:,cols); % swap indices based on solver generated function
    X = [X ones(length(X),1)]; % add extra column for material 
    LCOE = fval(:,1);
    Pvar = fval(:,2);
    lambdaActive=lambda.active';
    lambdaLower=lambda.lower';
    lambdaUpper=lambda.upper';

    new_objs = true; % switch between LCOE-Pvar and capex-Pavg
    
    [J1, minJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
     J2, minJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
     x_best_J1, x_best_J2, x_nom, x_balanced, idxo] = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w, new_objs);
    
    %% super simple "pareto" plot of just single objective optimizations
    showSingleObj = true;
    showImages = false;
    pareto_plot(J1, minJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, NaN*J1_balanced,...
                J2, minJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, NaN*J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, [], showSingleObj, showImages, p, new_objs)
    
    %% simple pareto plot
    showSingleObj = false;
    showImages = false;
    pareto_plot(J1, minJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
                J2, minJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj, showImages, p, new_objs)
    
    %% plot pareto front with annotations and embedded images of three recommended designs
    showSingleObj = true;
    showImages = true;
    pareto_plot(J1, minJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
                J2, minJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj, showImages, p, new_objs)
    
    %% plots for DVs as a fn of percent along the pareto
    J1_max = Inf;%p.LCOE_max;
    design_heuristics_plot(J1, minJ1, idx_best_J1, x_best_J1, ...
                           J2, minJ2, idx_best_J2, X, idxo, J1_max, b.var_names_pretty(1:end-1),new_objs)

    lagrange_multiplier_plot(lambdaActive,lambdaUpper,lambdaLower)

end
%%
function [J1, bestJ1, idx_best_J1, J1_nom, ...
         J1_nom_sim,  J1_solar,  J1_balanced,...
         J2, bestJ2, idx_best_J2, J2_nom, ...
         J2_nom_sim, J2_solar, J2_balanced,...
         x_best_J1, x_best_J2, x_nom, x_balanced, idxo] ...
                                        = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w, new_objs)
    
    if new_objs
        num_points = size(X,1);
        P_avg = zeros(num_points, 1);
        for i=1:num_points
            [~, ~, ~, ~, val] = simulation(X(i,:),p);
            P_avg(i) = val.power_avg;
        end
        J1 = P_avg;
        J1_fieldname = 'power_avg';
        J1_fieldidx = ':';
        J1_mult = 1/1000; % W to kW

        capex_design = Pvar; % second sim output means something different
        J2 = capex_design;
        J2_fieldname = 'capex_design';
        J2_fieldidx = 4;
        J2_mult = 1; % $M

        J1_solar = NaN;
        J2_solar = NaN;

        J2_balanced_approx = 300;
        min_max_sign = [-1 1]; % -1 for maximize, 1 for minimize

    else
        J1 = LCOE;
        J1_fieldname = 'LCOE';
        J1_fieldidx = 4;
        J1_mult = 1; % $/kWh

        J2 = Pvar;
        J2_fieldname = 'c_v';
        J2_fieldidx  = ':';
        J2_mult = 1; % percent (0-100)
    
        % solar
        J1_solar = 0.03; % LCOE solar
        J2_solar = 125;  % P var solar

        J2_balanced_approx = 100;
        min_max_sign = [1 1];
    end

    % find basic pareto front
    [~,idxo] = paretoFront(-min_max_sign .* [J1 J2]);
    [bestJ1,idx_best_J1] = min(min_max_sign(1)*J1);
    bestJ1 = min_max_sign(1) * bestJ1;
    [bestJ2,idx_best_J2] = min(min_max_sign(2)*J2);
    bestJ2 = min_max_sign(2) * bestJ2;

    % RM3 nominal - actual
    RM3_report = validation_inputs('report');
    RM3_wecsim = validation_inputs('wecsim');
    try
        J1_nom = [RM3_report.(J1_fieldname)(J1_fieldidx),...
                  RM3_wecsim.(J1_fieldname)(J1_fieldidx)];
        J2_nom = [RM3_report.(J2_fieldname)(J2_fieldidx),...
                  RM3_wecsim.(J2_fieldname)(J2_fieldidx)];
    catch
        J1_nom = [NaN NaN]; % fixme
        J2_nom = [NaN NaN];
    end
    
    % RM3 nominal - simulated
    x_nom = [b.X_noms b_w.X_noms];
    x_nom(end+1,:) = 1;

    [~, ~, ~, ~, val_report_sim] = simulation(x_nom(:,1),p);
    [~, ~, ~, ~, val_wecsim_sim] = simulation(x_nom(:,2),p_w);
    try
        J1_report_sim = val_report_sim.(J1_fieldname);
        J2_report_sim = val_report_sim.(J2_fieldname);
        J1_wecsim_sim = val_wecsim_sim.(J1_fieldname);
        J2_wecsim_sim = val_wecsim_sim.(J2_fieldname);
    catch
        J1_report_sim = NaN; % fixme
        J2_report_sim = NaN;
        J1_wecsim_sim = NaN;
        J2_wecsim_sim = NaN;
    end
    
    J1_nom_sim = [J1_report_sim, J1_wecsim_sim];
    J2_nom_sim = [J2_report_sim, J2_wecsim_sim];

    % balanced design
    [~,idx_balanced] = min(abs(J2-J2_balanced_approx));
    J1_balanced = J1(idx_balanced);
    J2_balanced = J2(idx_balanced);
    
    % idenitfy design variables for best designs
    x_best_J1 = X(idx_best_J1,:)
    x_best_J2 = X(idx_best_J2,:)
    x_balanced = X(idx_balanced,:)

    % apply multplier to all J1
    J1s = {J1, bestJ1, J1_nom, J1_nom_sim, J1_solar, J1_balanced};
    J1s_scaled = cellfun(@(x) x * J1_mult, J1s, 'UniformOutput', false);
    [J1, bestJ1, J1_nom, J1_nom_sim, ...
        J1_solar, J1_balanced] = deal(J1s_scaled{:});

    % apply multiplier to all J2
    J2s = {J2, bestJ2, J2_nom, J2_nom_sim, J2_solar, J2_balanced};
    J2s_scaled = cellfun(@(x) x * J2_mult, J2s, 'UniformOutput', false);
    [J2, bestJ2, J2_nom, J2_nom_sim, ...
        J2_solar, J2_balanced] = deal(J2s_scaled{:});
         
end

%%
function [] = pareto_plot(LCOE,minLCOE,idx_best_LCOE,LCOE_nom, LCOE_nom_sim, LCOE_solar, LCOE_balanced,...
                          Pvar,minPvar,idx_best_Pvar,P_var_nom,P_var_nom_sim,P_var_solar,P_var_balanced,...
                          x_best_LCOE,x_best_Pvar,x_nom,x_balanced,idxo,showSingleObj,showImages,p,new_objs)
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
    if new_objs
        xlabel('Average Electrical Power (kW)')
        ylabel('Structural and PTO Cost ($M)')
    else
        xlabel('LCOE ($/kWh)')
        ylabel('Power Variation (%)')
        xlim([0 1])
        ylim([25 165])
    end
    
    title('Pareto Front')
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
                                     overallX, idxo, LCOE_max, var_names, new_objs)

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
    
    if new_objs
        legend_text = {'Average Electrical Power (kW)','Structural and PTO Cost ($M)'};
        scale = 1;
    else
        scale = 100; % convert $ to cents
        cent = char(0162);
        legend_text = {['LCOE (' cent '/kWh)'],'c_v (%)'};
    end
    figure
    plot(pct_angle,scale*LCOE_pareto_sorted,'Color',[0.4940 0.1840 0.5560]) % purple
    hold on
    plot(pct_angle,Pvar_pareto_sorted,'Color',[0.4660 0.6740 0.1880]) % green
    grid on
    xlabel('Percent along the Pareto Curve')
    ylabel('Objective Value')
    improvePlot
    legend(legend_text,Location='northeast')
end

%%
function [] = lagrange_multiplier_plot(lambdaActive, ...
    lambdaUpper,lambdaLower);
figure
hold on
[m,n]=size(lambdaActive);
color_cell = {'.r','or','+r','*r','xr','-r','|r','+b','ob','+b','*b','xb','-b','|b'};
% I'll change these later
for i=1:m
    for j = 1:n
        color=color_cell{rem(j,15)};
        plot(i*100/60,lambdaActive(i,j),color,'MarkerSize',12);
    end
end
%plot(pct_angle,lambda_active_sorted)
title("Lagrange Multipliers")
xlabel("Percent Along the Pareto Front")
ylabel("Lagrange Multiplier")
legend("prevent float too heavy","prevent float too light", ...
    "prevent spar too heavy","prevent spar too light","stability",...
    "float survives max force","spar survives max force",...
    "damping plate survives max force","spar doesn't buckle", ...
    'positive power','damping plate diameter','prevent float rising above top of spar',...
    'prevent too expensive','prevent irrelevant max force')
hold off
%improvePlot

% lower bound plot
figure
hold on
[m,n]=size(lambdaLower);
for i=1:m
    for j = 1:n
        color=color_cell{rem(j,15)};
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
[m,n]=size(lambdaUpper);
for i=1:m
    for j = 1:n
        color=color_cell{rem(j,15)};
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
hold off
end