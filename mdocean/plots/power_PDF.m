function t = power_PDF(X,p,color,t,which_plots)

[~,~,P_matrix] = simulation(X, p);
P_matrix = P_matrix / 1e3; % convert W to kW
P_counts = repelem(P_matrix(:), round(p.JPD(:)*1000));

logscale = true; % toggle whether power is logscale
if logscale
    P_min = 10;
    P_max = 10000;
    prob_min = 1e-4;
    prob_max = 1;
    prob_label_pos = [10^-(1.5),1e-5];
else
    P_min = 50;
    P_max = 500;
    prob_min = 0;
    prob_max = 0.005;
    prob_label_pos = [-5 -.0015];
end

if nargin < 5
    which_plots = [true true];
end

if which_plots(1) % PDF
    ax1 = nexttile(1);
    fit = true; % false for histogram, true for ksdensity fit
    if ~fit
%         h = histogram(P_counts,'BinWidth',10,'Normalization','pdf');
        
        [N,edges] = histcounts(P_counts,'BinWidth',.1,'Normalization','pdf');
        edge_centers = 1/2 *(edges(1:end-1) + edges(2:end));
        edges = 10.^(-1:0.1:3); 
        hist = histc(P_counts, edges);
        centers = sqrt(edges(1:end-1).*edges(2:end));
        bar(edges,hist)
        h = plot(edge_centers,N,'Color',color);
        hold on
    else
    
         if logscale
             P_x = logspace(log10(P_min),log10(P_max),100);
         else
             P_x = linspace(P_min,P_max,100);
         end
        [f,xi]= ksdensity(P_counts,P_x,'support', 'positive');
        if logscale
            h = semilogx(xi,f);
            ax1.XScale = 'log';
            ax1.YScale = 'log';
            ax1.YMinorGrid = 'off';
        else
            h = plot(xi,f);
        end
    end
    ylim([prob_min,prob_max])
        
    grid on
    if nargin > 2
        h.Color = color;
    end
    hold on
    title('Probability Density Function')
    yy = ylabel('Probability (-)');
    yy.Position = [prob_label_pos 0];
end

if which_plots(2) % CDF
    ax2 = nexttile(2);
    h = cdfplot(P_counts);
    hold on
    if nargin > 2
        h.Color = color;
    end
    title('Cumulative Distribution Function')
    ylabel('Portion of Hours per Year')
    
    if logscale
        ax2.XScale = 'log';
    end
    grid on
    
    xlabel('Average Electrical Power (kW)')
end

% align axes if plotting both
if which_plots(1) && which_plots(2)
    linkaxes([ax1,ax2],'x');
    xticklabels(ax1,{})
    xticks(ax1,xticks(ax2))
end

xlim([P_min P_max])
improvePlot

nexttile(1)


end