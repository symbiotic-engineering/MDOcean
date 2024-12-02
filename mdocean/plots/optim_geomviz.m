function flag = optim_geomviz(x,p,b)

    X =[x(b.idxs_recover); 1]; % reorder indices

    grey = [0.8, 0.8, 0.8]; % RGB for light grey

    % Get all line and rectangle objects from previous iterations
    lines = findobj(gca, 'Type', 'line');
    rectangles = findobj(gca, 'Type', 'rectangle');
    
    % Change their color to light grey
    for i = 1:length(lines)
        lines(i).Color = grey;
    end
    for i = 1:length(rectangles)
        rectangles(i).EdgeColor = grey;
    end

    % plot geomtry of current iteration
    visualize_geometry(X,p,true)

    flag = false; % don't stop optimization

end