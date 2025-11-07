function stop = optimplotx_custom(x, optimValues, state, b)
    % optimplotfvalconstr Plot value of the objective function vs iteration.
    % Infeasible points are marked red. If no fval is available, this function
    % will plot constraint violation (infeasibility).
    %

    %   Copyright 2019-2020 The MathWorks, Inc.

    persistent plotBest feasLegend

    stop = false;

    if strcmpi(state, 'init')
        plotBest = gobjects(1, length(x));
        feasLegend = [];
    end

    if optimValues.funccount == 0 || (isempty(optimValues.fval) && isempty(optimValues.constrviolation))
        % no function evals or none of the trials are successfully evaluated; no plots.
        return
    end

    X = [x(b.idxs_recover); 1]; % reorder indices
    best = X;

    marks = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

    for i = 1:length(x)
        if isa(plotBest(i), 'matlab.graphics.GraphicsPlaceholder')
            xlabel(getString(message('MATLAB:optimfun:funfun:optimplots:LabelIteration')), 'interp', 'none');

            if isempty(optimValues.fval)
                ylabel_msg = getString(message('optim:optimplot:YlabelConstrViol'));
                title_msg = getString(message('optim:optimplot:TitleMaxConstrViol', sprintf('%g', best)));
            else
                ylabel_msg = getString(message('MATLAB:optimfun:funfun:optimplots:LabelFunctionValue'));
                title_msg = getString(message('optim:optimplot:TitleBestFunVal', sprintf('%g', best)));
            end
            ylabel(ylabel_msg, 'interp', 'none');
            title(title_msg, 'interp', 'none');
            hold on;
            grid on;
        end

        % Plot points
        if isa(plotBest(i), 'matlab.graphics.GraphicsPlaceholder')
            plotBest(i) = plot(optimValues.iteration, best(i), marks{i}, 'DisplayName', b.var_names_pretty{i});
            set(plotBest(i), 'Tag', 'plotbestf');
            %         legendHndl(end+1) = plotBest;
            %         legendStr{end+1} = getString(message('optim:optimplot:LegendBestFval'));
        else
            newX = [get(plotBest(i), 'Xdata') optimValues.iteration];
            newY = [get(plotBest(i), 'Ydata') best(i)];
            set(plotBest(i), 'Xdata', newX, 'Ydata', newY);
        end

    end

    title_msg = 'Current x';
    set(get(gca, 'Title'), 'String', title_msg, 'interp', 'none');

    if isempty(feasLegend)
        feasLegend = legend();
    end

    if strcmp(state, 'done')
        hold off;
    end
end
