function f = color_each_element(matrix)

f = figure;
imagesc(matrix);
colormap(bluewhitered);

ax = gca;

nx = size(matrix,2);
ny = size(matrix,1);

set(ax,'XTick',-.5 + 1:nx)
set(ax,'YTick',-.5 + 1:ny)
set(ax,'XTickLabel',1:nx)
set(ax,'YTickLabel',1:ny)

grid on;
axis image

colorbar;

end
