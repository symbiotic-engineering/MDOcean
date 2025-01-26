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

    new_objs = true; % switch between LCOE-Pvar and capex-Pavg
    
    [J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
     J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
     x_best_J1, x_best_J2, x_nom, x_balanced, idxo] = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w, new_objs);
    
    %% super simple "pareto" plot of just single objective optimizations
    showSingleObj = true;
    showImages = false;
    pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, NaN*J1_balanced,...
                J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, NaN*J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, [], showSingleObj, showImages, p, new_objs)
    
    %% simple pareto plot
    showSingleObj = false;
    showImages = false;
    pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim*NaN, J1_solar, J1_balanced,...
                J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim*NaN, J2_solar, J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj, showImages, p, new_objs)
    
    %% plot pareto front with annotations and embedded images of three recommended designs
    showSingleObj = true;
    showImages = true;
    pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
                J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj, showImages, p, new_objs)
    
    %% plots for DVs as a fn of percent along the pareto
    J1_max = Inf;%p.LCOE_max;
    design_heuristics_plot(J1, bestJ1, idx_best_J1, x_best_J1, ...
                           J2, bestJ2, idx_best_J2, X, idxo, J1_max, b.var_names_pretty(1:end-1),new_objs)

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
function [] = pareto_plot(J1,bestJ1,idx_best_J1,J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
                          J2,bestJ2,idx_best_J2,J2_nom,J2_nom_sim,J2_solar,J2_balanced,...
                          x_best_J1,x_best_J2,x_nom,x_balanced,idxo,showSingleObj,showImages,p,new_objs)
    figure
    % overall pareto front
    plot(J1(idxo),J2(idxo),'bs','MarkerFaceColor','b','HandleVisibility','off')
    hold on
    
    % utopia point
    plot(bestJ1,bestJ2,'gp','MarkerFaceColor','g','MarkerSize',20,'HandleVisibility','off')
    
    % RM3 nominal reference - report
    plot(J1_nom(1),J2_nom(1),'rd','HandleVisibility','off')
    plot(J1_nom_sim(1),J2_nom_sim(1),'rs','HandleVisibility','off')

    % RM3 nominal reference - wecsim
    plot(J1_nom(2),J2_nom(2),'md','HandleVisibility','off')
    plot(J1_nom_sim(2),J2_nom_sim(2),'ms','HandleVisibility','off')
    
    if showSingleObj  
        % black squares for 3 ref points
        plot(bestJ1,J2(idx_best_J1),'ks','HandleVisibility','off')
        plot(J1(idx_best_J2),bestJ2,'ks','HandleVisibility','off')
        plot(J1_balanced,J2_balanced,'ks','HandleVisibility','off')
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
    plot(J1_solar, J2_solar,'o','MarkerSize',12,'MarkerEdgeColor',yellow,...
        'MarkerFaceColor',yellow,'HandleVisibility','off')
    % for the yellow color to work, do not use improvePlot below here
    
    % text labels
    sz = 14;
    text(J1_nom(1)+10,J2_nom(1),'RM3 Report','FontSize',sz)
    text(J1_nom(1)+.01,J2_nom(1)-5,'Actual [10]','FontSize',sz)
    text(J1_nom_sim(1)-.02,J2_nom_sim(1)+5,'RM3 Report','FontSize',sz)
    text(J1_nom_sim(1)+.03,J2_nom_sim(1),'Sim','FontSize',sz)

    text(J1_nom(2)+.03,J2_nom(2),'RM3 WecSim','FontSize',sz)
    text(J1_nom(2)+.01,J2_nom(2)-5,'Actual','FontSize',sz)
    text(J1_nom_sim(2)-.02,J2_nom_sim(2)+5,'RM3 WecSim','FontSize',sz)
    text(J1_nom_sim(2)+.03,J2_nom_sim(2),'Sim','FontSize',sz)

    text(bestJ1+.03,bestJ2-2,'Utopia Point','FontSize',sz)
    text(J1_solar+.03,J2_solar,'Solar','FontSize',sz)
    if showSingleObj
        text(J1(idx_best_J1)+.03,J2(idx_best_J1),'Cheapest','FontSize',sz)
        text(J1(idx_best_J2)+.03,J2(idx_best_J2)-3,'Least Variable','FontSize',sz)
        text(J1_balanced-.15,J2_balanced+5,'Balanced Design','FontSize',sz)
    end

    showLCOEContours = true;
    if showLCOEContours
        LCOE_min = 0.083;
        LCOE_nom = 0.76; % fixme hardcoded
        overlay_LCOE(p, LCOE_nom, LCOE_min)
    end

    if showImages   % small pictures of best geometries

        upper_left = [.28 .6];
        mini_plot(upper_left,x_best_J1,p)

        lower_right = [.51 .23];
        mini_plot(lower_right,x_best_J2,p)
        
        balanced_pos = [.10 .28];
        mini_plot(balanced_pos,x_balanced,p)

        report_pos = [.7 .53];
        mini_plot(report_pos,x_nom(:,1),p)
        
        wecsim_pos = [.8 .7];
        mini_plot(wecsim_pos,x_nom(:,2),p)
    end

end

function mini_plot(pos, x, p)
    mini_plot_size = [.2 .22];
    axes('Position',[pos mini_plot_size])
    box on
    visualize_geometry(x,p,true);
    set(gca,'XTickLabel',[],'YTickLabel',[])
end

%%
function [] = design_heuristics_plot(overallJ1, bestJ1, idx_best_J1, x_best_J1, ...
                                     overallJ2, bestJ2, idx_best_J2, ...
                                     overallX, idxo, J1_max, var_names, new_objs)

    J1_pareto = overallJ1(idxo);
    J2_pareto = overallJ2(idxo);
    [J1_pareto_sorted,idx_sort] = sort(J1_pareto(J1_pareto<J1_max));
    J2_pareto_sorted = J2_pareto(J1_pareto<J1_max);
    J2_pareto_sorted = J2_pareto_sorted(idx_sort);
    
    pct = linspace(0,100,length(idx_sort));

    J1_worst = overallJ1(idx_best_J2);
    J2_worst = overallJ2(idx_best_J1);
    J1_range = bestJ1 - J1_worst;
    J2_range = bestJ2 - J2_worst;
    frac_J1 = abs( (J1_worst - J1_pareto_sorted) / J1_range );
    frac_J2 = abs( (J2_worst - J2_pareto_sorted) / J2_range );

    % this keeps 0% corresponding to left (min J1) and 100% corresponding to right (max J1)
    if new_objs
        dir = {'max','min'};
    else
        dir = {'min','min'};
    end
    if strcmp(dir(1),'max')
        num = frac_J1;
        den = frac_J2;
    elseif strcmp( dir(1),'min')
        num = frac_J2;
        den = frac_J1;
    else
        error('invalid direction string')
    end

    pct_angle = 100/(pi/2) * atan( num ./ den);
    pct_angle(pct_angle==-100) = 100;
    
    X_pareto = overallX(idxo,:);
    X_pareto = X_pareto(J1_pareto<J1_max,:);
    X_pareto_sorted = X_pareto(idx_sort,:);
    
    X_pareto_sorted_scaled = X_pareto_sorted ./ repmat(x_best_J1,length(idx_sort),1);
    
    X_pareto_sorted_scaled = X_pareto_sorted_scaled(:,1:end-1); % get rid of material
    
    windowSize = round(length(idx_sort) * 5/100);
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    n_points_interp = 51;
    pct_angle_even = linspace(0,100,n_points_interp);
    y = zeros(n_points_interp, size(X_pareto_sorted_scaled,2));
    
    for i=1:size(X_pareto_sorted_scaled,2)
        x = X_pareto_sorted_scaled(:,i);
        x_even_spread = interp1(pct_angle,x,pct_angle_even).';
        if strcmp(dir(1),'max')
            x_even_spread = flipud(x_even_spread); % apply filter backwards
        end
        x_padded = [ones(windowSize,1); x_even_spread];
        yy = filter(b,a,x_padded); % moving average filter
        yy_crop = yy((windowSize+1):end);
        if strcmp(dir(1),'max')
            yy_crop = flipud(yy_crop); % apply filter backwards
        end
        y(:,i) = yy_crop;
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
    cols = {'r:','r--','r-','r-.','r.',...       % bulk dims
            'b:','b--',...                       % PTO
            'g:','g--','g-','g-.','g.'}; % 'g*'  % structural
    figure
    for i=1:size(X_pareto_sorted_scaled,2)
        semilogy(pct_angle_even,y(:,i),cols{i})
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
    ylim([.03 2])
    set(gca,'YTick',y_tick)
    improvePlot
    legend(var_names,'Location','eastoutside')
    set(gca,'YMinorGrid','on')
    set(gca,'XGrid','on')
    set(gca, 'Children', flipud(get(gca, 'Children')) ) % put fake gridlines behind real lines
    
    if new_objs
        legend_text = {'Average Electrical Power','Structural and PTO Cost'};
        scale_1 = 1;
        scale_2 = 1; %100;
    else
        scale_1 = 100; % convert $ to cents
        scale_2 = 1;
        cent = char(0162);
        legend_text = {['LCOE (' cent '/kWh)'],'c_v (%)'};
    end
    figure
    plot(pct_angle,scale_1*J1_pareto_sorted/bestJ1,'Color',[0.4940 0.1840 0.5560]) % purple
    hold on
    plot(pct_angle,scale_2*J2_pareto_sorted/J2_worst,'Color',[0.4660 0.6740 0.1880]) % green
    grid on
    xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Objective Value')
    improvePlot
    legend(legend_text,Location='northeast')
end

function [] = overlay_LCOE(p, LCOE_nom, LCOE_min)
    n_points = 30;
    P_elec_kW = xlim;
    C_design_M = ylim;
    P_elec_W = 1000*linspace(P_elec_kW(1), P_elec_kW(2), n_points);
    C_design = 1e6*linspace(C_design_M(1), C_design_M(2), n_points);

    [P_ELEC,C_DESIGN] = meshgrid(P_elec_W,C_design);
    LCOE = LCOE_fcn(P_ELEC, C_DESIGN, p);

    hold on
    levs = [.06;
            LCOE_min;
            .1; .15; .2; .3; .5; .8];
            %LCOE_nom];
    grey = [.85 .85 .85];
    [c,h] = contour(P_ELEC/1000,C_DESIGN/1e6,LCOE,levs,'Color',grey);
    clabel(c,h,'manual');
    legend('LCOE')

end

function LCOE = LCOE_fcn(P_elec, C_design, p)
% fixme I've just hardcoded stuff from econ, this could become inconsistent easily
% this fcn is vectorized, whereas econ is not.
    development     = 4553000;
    infrastructure  = 990000;
    mooring         = p.N_WEC * 525000;
    profitmargin    = 356000;
    installation    = 5909000;
    contingency     = 1590000;
    
    capex = development + infrastructure + mooring + C_design * p.N_WEC ...
            + profitmargin + installation + contingency; 
    
    % opex reflects cost of operation, postinstall, replacement, consumables, and
    % insurance
    % opex equation was generated by Casio online power regression calculator
    opex = 1193000 * p.N_WEC^0.4433;
    
    hr_per_yr = 8766;
    P_avg = p.N_WEC * P_elec * p.eff_array;
    aep = P_avg * hr_per_yr / 1000; % annual energy production: W to kWh per year
    
    LCOE = (p.FCR*capex + opex)./aep;
end