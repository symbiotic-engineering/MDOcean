add_mdocean_path;
clear all;
mathfonts = {'Noto Sans Math','STIX Math','Latin Modern Math','TeX Gyre Termes Math','TeX Gyre Bonum Math', 'DejaVu Math TeX Gyre','TeX Gyre Pagella Math','TeX Gyre Schola Math','TeX Gyre DejaVu Math','DejaVu Sans','Noto Sans','Noto Sans Symbols'};

% 'Asana Math',

fonts = [listfonts; mathfonts.'];

fig=figure;

% must set limits before raster_text to get correct scaling
xlim([0   length(fonts)+1]);
ylim([0 1]);

for i=1:length(fonts)
    text(i, .5, unicode('1d54f'),'FontName',fonts{i},'Interpreter','none'); 
    hold on
    raster_text(gca, i, .1, [unicode('1d54f') '_L'], .5/length(fonts), fonts{i});
end

exportgraphics(fig,'test_font.pdf');