function saveSVG(filename,svg)


fid=fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n');
fprintf(fid,'<!-- Created with Matlab  -->\n');
fprintf(fid,'<svg\n');
fprintf(fid,'   xmlns:dc="http://purl.org/dc/elements/1.1/"\n');
fprintf(fid,'   xmlns:cc="http://creativecommons.org/ns#"\n');
fprintf(fid,'   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n');
fprintf(fid,'   xmlns:svg="http://www.w3.org/2000/svg"\n');
fprintf(fid,'   xmlns="http://www.w3.org/2000/svg"\n');
fprintf(fid,'   xmlns:xlink="http://www.w3.org/1999/xlink"\n');
fprintf(fid,'   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"\n');
fprintf(fid,'   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n');
fprintf(fid,'   width="2007.2064"\n');
fprintf(fid,'   height="3009.4392"\n');
fprintf(fid,'   id="svg2"\n');
fprintf(fid,'   version="1.1">\n');

fprintf(fid,' <g\n');
fprintf(fid,'     inkscape:groupmode="layer"\n');
fprintf(fid,'     id="layer2"\n');
fprintf(fid,'     inkscape:label="image"\n');
fprintf(fid,'     style="display:inline">\n');

for  k=1:length(svg.images)
    image=svg.images{k};
    fprintf(fid,'   <image\n');
    fprintf(fid,'      sodipodi:absref="%s"\n',image.file);
    fprintf(fid,'      xlink:href="building.jpg"\n');
    fprintf(fid,'      y="%d"\n',image.y);
    fprintf(fid,'      x="%d"\n',image.x);
    fprintf(fid,'      id="image3063"\n');
    fprintf(fid,'      height="%d"\n',image.height);
    fprintf(fid,'      width="%d" />\n',image.width);
    
end
fprintf(fid,'</g>\n');

for  k=1:length(svg.layers)
    layer=svg.layers{k};
    nbpolys=length(layer.polys);
    
    fprintf(fid,'<g\n');
    fprintf(fid,'  inkscape:groupmode="layer"\n');
    fprintf(fid,'  id="layer%d"',k);
    fprintf(fid,'  inkscape:label="%s"\n', layer.name);
    fprintf(fid,'  style="display:inline">\n');
    
    for idpoly=1:nbpolys
        
        fprintf(fid,'<path ');
        colorhex=color2hex(layer.colors(:,idpoly));
        strokecolorhex=color2hex(layer.stroke_colors(:,idpoly));
        fprintf(fid,'   style="opacity:%f;;fill:#%s;fill-opacity:1;stroke:#%s;stroke-width:1;stroke-miterlimit:4;stroke-opacity:1;stroke-dasharray:none" \n',layer.colors(4,idpoly),colorhex,strokecolorhex);
       d=[layer.polys{idpoly}(:,1),diff(layer.polys{idpoly},1,2)];
        s=sprintf('%f,%f  ', d);
        fprintf(fid,'   d="m %s z"\n',s);
        if   not(isempty(layer.svgids{idpoly}))
            fprintf(fid,'   id="%s"\n',layer.svgids{idpoly});
        end
        fprintf(fid,'   inkscape:connector-curvature="0"\n');
        fprintf(fid,'   sodipodi:nodetypes="ccccc" />\n');  
    end
    fprintf(fid,'</g>\n');
end

fprintf(fid,'</svg>\n');
fclose(fid);
function hex=color2hex(color)
h1=dec2hex(round(color(1)*255),2);
h2=dec2hex(round(color(2)*255),2);
h3=dec2hex(round(color(3)*255),2);
hex=lower(sprintf( '%s%s%s',h1,h2,h3));

