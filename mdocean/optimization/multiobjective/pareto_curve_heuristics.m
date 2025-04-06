function pareto_curve_heuristics()
    p0 = parameters();
    b = var_bounds();
    p_w = parameters('wecsim');
    b_w = var_bounds('wecsim');
    
    d=dir("**/pareto_search_results*");
    load(d(end).name)

    if ~isequaln(p,p0)
        warning(['You are loading results with different parameters than your ' ...
            'local machine right now. WecSim validation results (p_w) may be incorrect.'])
    end

    if ~exist('tol','var')
        tol = 1e-6;
    end

    new_objs = true; % switch between LCOE-Pvar and capex-Pavg

    constraint_active_plot(residuals,fval,tol,b,new_objs);

    cols = b.idxs_recover;
    X = x(:,cols); % swap indices based on solver generated function
    X = [X ones(length(X),1)]; % add extra column for material 
    LCOE = fval(:,1);
    Pvar = fval(:,2);

    [J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
     J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
     x_best_J1, x_best_J2, x_nom, x_balanced, idxo, LCOE_nom] = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w,new_objs);
    
    %% super simple "pareto" plot of just single objective optimizations
    showSingleObj = true;
    showImages = false;
    showLCOEContours = false;
    pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, NaN*J1_balanced,...
                J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, NaN*J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, [], showSingleObj, ...
                showImages, showLCOEContours, p, new_objs)
    
    %% simple pareto plot
    showSingleObj = false;
    showImages = false;
    showLCOEContours = true;
    pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim.*[1 NaN], J1_solar, J1_balanced,...
                J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim.*[1 NaN], J2_solar, J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj, ...
                showImages, showLCOEContours, p, new_objs, LCOE_nom, min(LCOE))
    
    %% plot pareto front with annotations and embedded images of three recommended designs
    showSingleObj = true;
    showImages = true;
    showLCOEContours = false;
    pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim.*[1 NaN], J1_solar, J1_balanced,...
                J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim.*[1 NaN], J2_solar, J2_balanced,...
                x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj, ...
                showImages, showLCOEContours, p, new_objs)
    
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
         x_best_J1, x_best_J2, x_nom, x_balanced, ...
         idxo, LCOE_nom] = process_pareto_front(LCOE,Pvar,X,p,p_w,b,b_w,new_objs)
    
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
        J2_fieldname = 'J_capex_design';
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
        LCOE_nom = RM3_report.LCOE(4);
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
function [] = pareto_plot(J1, bestJ1, idx_best_J1, J1_nom, J1_nom_sim, J1_solar, J1_balanced,...
                          J2, bestJ2, idx_best_J2, J2_nom, J2_nom_sim, J2_solar, J2_balanced,...
                          x_best_J1, x_best_J2, x_nom, x_balanced, idxo, showSingleObj,...
                          showImages, showLCOEContours, p, new_objs, LCOE_nom, LCOE_min)
    figure
    % overall pareto front
    plot(J1(idxo),J2(idxo),'bs','MarkerFaceColor','b','HandleVisibility','off')
    hold on
    
    % utopia point
    plot(bestJ1,bestJ2,'gp','MarkerFaceColor','g','MarkerSize',20,'HandleVisibility','off')
    
    % RM3 nominal reference - report
    plot(J1_nom(1),    J2_nom(1),    'rd','HandleVisibility','off')
    plot(J1_nom_sim(1),J2_nom_sim(1),'rs','HandleVisibility','off')

    % RM3 nominal reference - wecsim
    plot(J1_nom(2),    J2_nom(2),    'md','HandleVisibility','off')
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
        xlim([.9*min(J1) 1.1*max(J1)])
        ylim([.9*min(J2) 1.05*max(J2)])
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
    text(J1_nom(1)+2,J2_nom(1)+.02,'RM3 Report','FontSize',sz)
    text(J1_nom(1)+.01,J2_nom(1)-5,'Actual [10]','FontSize',sz)
    text(J1_nom_sim(1)-.02,J2_nom_sim(1)+5,'RM3 Report','FontSize',sz)
    text(J1_nom_sim(1)+.03,J2_nom_sim(1)-.05,'Sim','FontSize',sz)

    text(J1_nom(2)+.03,J2_nom(2),'RM3 WecSim','FontSize',sz)
    text(J1_nom(2)+.01,J2_nom(2)-5,'Actual','FontSize',sz)
    text(J1_nom_sim(2)-.02,J2_nom_sim(2)+5,'RM3 WecSim','FontSize',sz)
    text(J1_nom_sim(2)+.03,J2_nom_sim(2),'Sim','FontSize',sz)

    text(bestJ1+.03,bestJ2-2,'Utopia Point','FontSize',sz)
    text(J1_solar+.03,J2_solar,'Solar','FontSize',sz)
    if showSingleObj
        text(J1(idx_best_J1)+.03,J2(idx_best_J1),'Most Power','FontSize',sz)
        text(J1(idx_best_J2)+.03,J2(idx_best_J2)-3,'Least Variable','FontSize',sz)
        text(J1_balanced-.15,J2_balanced+5,'Balanced Design','FontSize',sz)
    end

    if showLCOEContours
        overlay_LCOE(p, LCOE_nom, LCOE_min)
    end

    if showImages   % small pictures of best geometries
        ylim([0.64, 2.55])

        upper_left = [.75 .63]; %.28,.6
        mini_plot(upper_left,x_best_J1,p)

        lower_right = [.09 .19]; %.51,.23
        mini_plot(lower_right,x_best_J2,p)
        
        balanced_pos = [.42 .12]; %0.1,.35
        mini_plot(balanced_pos,x_balanced,p)

        report_pos = [.09 .55]; %.7,.53
        mini_plot(report_pos,x_nom(:,1),p)
        
        %wecsim_pos = [.28 .65]; %.8,.7
        %mini_plot(wecsim_pos,x_nom(:,2),p)
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
function [] = design_heuristics_plot(overallJ1, J1_best, idx_best_J1, x_best_J1, ...
                                     overallJ2, J2_best, idx_best_J2, ...
                                     overallX, idxo, J1_max, var_names, new_objs)

    % filter out pareto nonoptimal points
    J1_pareto = overallJ1(idxo);
    J2_pareto = overallJ2(idxo);
    X_pareto  = overallX(idxo,1:end-1); % get rid of material

    % filter out J1 above thresh
    idx_keep = J1_pareto < J1_max;

    J1_worst = overallJ1(idx_best_J2);
    J2_worst = overallJ2(idx_best_J1);

    % desired direction
    if new_objs
        dir = {'max','min'};
    else
        dir = {'min','min'};
    end

    % order from left to right on pareto front appearance, 
    % and normalize by optimal J1 values
    [X_pareto_ordered_normed,...
     J1_pareto_ordered_normed, ...
     J2_pareto_ordered_normed,...
     pct_angle] = pareto_order_normalize(J1_pareto(idx_keep), J2_pareto(idx_keep), ...
                                         X_pareto(idx_keep,:), dir,...
                                         J1_best, J1_worst, ...
                                         J2_best, J2_worst, x_best_J1(1:end-1));
    
    % apply filter
    backwards = strcmp(dir(1),'max');
    [X_filtered,pct_angle_even] = low_pass_filter(pct_angle, X_pareto_ordered_normed, backwards);
    
    % plot unfiltered
    figure
    semilogy(pct_angle,X_pareto_ordered_normed)
    title('Unfiltered Design Heuristics')
    xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Optimal Design Value')
    legend(var_names,'Location','eastoutside')
    ylim([0 15])
    improvePlot
    grid on
    set(gca,'YMinorGrid','on')
    
    % plot filtered
    cols = {'r:','r--','r-','r-.','r.',...       % bulk dims
            'b:','b--',...                       % PTO
            'g:','g--','g-','g-.','g.'}; % 'g*'  % structural
    figure
    s1 = subplot(2,1,1);
    for i=1:size(X_pareto_ordered_normed,2)
        semilogy(pct_angle_even,X_filtered(:,i),cols{i})
        hold on
    end
    
    % make fake major grid lines (didn't do "grid on" because then major lines show up for .2 .5 2 5 too)
    x_grid = [3 97];
    y_grid = [.1 1];
    y_tick = [.01 .02 .05 .1 .2 .5 1 2 5];
    for i=1:length(y_grid)
        plot(x_grid, y_grid(i)*[1 1],'Color',[.85 .85 .85]);
    end
    
    title('Design Trends')
    %xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Optimal Design Value')
    ylim([.9*min(X_filtered,[],'all'), 1.1*max(X_filtered,[],'all')])
    set(s1,'YTick',y_tick)
    improvePlot
    
    l1 = legend(var_names,'Location','eastoutside');
    set(s1,'YMinorGrid','on')
    set(s1,'XGrid','on')
    set(s1, 'Children', flipud(get(s1, 'Children')) ) % put fake gridlines behind real lines
    
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

    % objective heuristics
    s2 = subplot(2,1,2); %figure
    plot(pct_angle, scale_1 * J1_pareto_ordered_normed,'Color',[0.4940 0.1840 0.5560]) % purple
    hold on
    plot(pct_angle, scale_2 * J2_pareto_ordered_normed,'Color',[0.4660 0.6740 0.1880]) % green
    grid on
    xlabel('Percent along the Pareto Curve')
    ylabel('Normalized Objective Value')
    improvePlot
    l2 = legend(legend_text,Location='eastoutside');

    set(gcf,"Position",[100 100 715 835]) % make wider so x major grid shows for every 20
    l1.Position(1:2) = [.78 .43];
    s1.Position(2:4) = [.5 .6 .4];
    s2.Position(3:4) = [.6 .28];
    l2.Position(1:2) = [.54 .16];
end

function [filtered, t_even_spread] = low_pass_filter(time, signals, backwards)

    num_sigs = size(signals,2);

    windowSize = round(length(time) * 5/100);
    b = (1/windowSize) * ones(1,windowSize);
    a = 1;

    n_points_interp = 51;
    t_even_spread = linspace(min(time),max(time),n_points_interp);

    filtered = zeros(n_points_interp, num_sigs);
    
    for i=1:num_sigs
        % interpolate signal onto uniform time grid
        sig = signals(:,i);
        sig_even_spread = interp1(time,sig,t_even_spread).';

        % flip if filtering backwards in time
        if backwards
            sig_even_spread = flipud(sig_even_spread);
        end

        % pad, apply moving average filter, and crop away padding
        sig_padded = [ones(windowSize,1); sig_even_spread];
        yy = filter(b,a,sig_padded); 
        yy_crop = yy((windowSize+1):end);

        % flip if filtering backwards in time
        if backwards
            yy_crop = flipud(yy_crop); % apply filter backwards
        end

        filtered(:,i) = yy_crop;
    end
end

function [] = overlay_LCOE(p, LCOE_nom, LCOE_min)
    n_points = 30;
    P_elec_kW = xlim;
    C_design_M = ylim;
    P_elec_W = 1000*linspace(P_elec_kW(1), P_elec_kW(2), n_points);
    C_design = 1e6*linspace(C_design_M(1), C_design_M(2), n_points);

    [P_ELEC,C_DESIGN] = meshgrid(P_elec_W,C_design);
    LCOE = LCOE_from_capex_design_power(C_DESIGN, p.N_WEC, P_ELEC, p.FCR, p.eff_array);

    hold on
    levs = [LCOE_min;
            .12; .15; .2; .3; .5;
            LCOE_nom];
    grey = [.85 .85 .85];
    [~,h] = contour(P_ELEC/1000,C_DESIGN/1e6,LCOE,levs,'Color',grey,'ShowText','on','LabelSpacing',400);
    if isMATLABReleaseOlderThan("R2022b")
        h.LevelList = round(h.LevelList,3); % for old versions, round the contours directly. Will be slightly off from scatter points.
    else
        h.LabelFormat = '%0.2f'; % for new versions, just round the text
    end

    legend('LCOE ($/kWh)','Location','northwest')

end

function [X_pareto_ordered_normed,...
          J1_pareto_ordered_normed, ...
          J2_pareto_ordered_normed, pct_angle] = pareto_order_normalize(J1_pareto, J2_pareto, X_pareto, dir,  ...
                                                                        J1_best, J1_worst, J2_best, J2_worst,...
                                                                        X_best_J1)

    % sort from min J1 to max J1 (left to right)
    [J1_pareto_ordered,idx_order] = sort(J1_pareto);
    J2_pareto_ordered = J2_pareto(idx_order);
    X_pareto_ordered = X_pareto(idx_order,:);
    
    % normalize by best J1
    X_pareto_ordered_normed = X_pareto_ordered ./ repmat(X_best_J1,length(idx_order),1);
    J1_pareto_ordered_normed = J1_pareto_ordered / J1_best;
    J2_pareto_ordered_normed = J2_pareto_ordered / J2_worst;

    % ranges
    J1_range = J1_best - J1_worst;
    J2_range = J2_best - J2_worst;

    % fraction along each objective
    frac_J1 = abs( (J1_worst - J1_pareto_ordered) / J1_range );
    frac_J2 = abs( (J2_worst - J2_pareto_ordered) / J2_range );

    % choose num and den to keep 0% corresponding to left (min J1) and 
    % 100% corresponding to right (max J1). See notebook p83 1/25/25
    if strcmp(dir(1),'max')
        num = frac_J1;
        den = frac_J2;
    elseif strcmp( dir(1),'min')
        num = frac_J2;
        den = frac_J1;
    else
        error('invalid direction string')
    end

    % find angle along circle with origin at nadir point - useful because 
    % diminishing returns often makes pareto look like a quarter circle
    angle = atan( num ./ den);

    % express angle as percentage: 0 to pi/2 maps 0 to 100
    pct_angle = 100/(pi/2) * angle;
    pct_angle(pct_angle==-100) = 100;
    
end
