function h = raster_text(ax, x, y, str, x_span, fontname)
% creates a rasterized version of the specified text and places it on the given axes at (x,y) with the specified font. 
% x_span controls the width of the text as a fraction of the axis width. 
% This is useful for unicode characters that may not render correctly in vector format.
% You should set the axis limits before calling this function to ensure correct scaling.

    % Create offscreen figure and render glyph
    f = figure('Visible','off','Color','white');
    axf = axes(f,'Position',[0 0 1 1],'Visible','off');
    text(0.5,0.5,str,'FontSize',250,'FontName',fontname,...
        'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex');
    axis(axf,'off')
    drawnow

    % Export and read PNG
    imgFile = 'tmp_symbol.png';
    exportgraphics(f, imgFile, 'BackgroundColor','white');
    close(f)

    img = imread(imgFile);

    % Get axis data ranges and pixel sizes
    xLim = xlim(ax);
    yLim = ylim(ax);
    xRange = diff(xLim);
    yRange = diff(yLim);
    axPos = getpixelposition(ax);
    axPixelWidth = axPos(3);
    axPixelHeight = axPos(4);

    % Image dimensions and aspect ratio
    [imgH, imgW] = size(img, [1 2]);
    imgAspect = imgH / imgW;

    % Desired data-space width (20% of axis width)
    desiredDataWidth = x_span*xRange;

    % Convert to pixels and compute corresponding height
    px_per_x = axPixelWidth / xRange;
    desiredPixelsWidth = desiredDataWidth * px_per_x;
    desiredPixelsHeight = desiredPixelsWidth * imgAspect;

    % Convert back to data space
    px_per_y = axPixelHeight / yRange;
    desiredDataHeight = desiredPixelsHeight / px_per_y;

    xdata = [x - desiredDataWidth/2, x + desiredDataWidth/2];
    ydata = [y + desiredDataHeight/2, y - desiredDataHeight/2];

    % Place image on caller axes
    holdState = ishold(ax);
    hold(ax,'on')
    h = image(ax, xdata, ydata, img, 'Clipping','off');
    set(ax,'YDir','normal')
    if ~holdState
        hold(ax,'off')
    end

    delete(imgFile);
end