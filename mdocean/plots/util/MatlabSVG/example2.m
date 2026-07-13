svg=loadSVG('buildingNoImage.svg');
svg;
% svg.images{1}
% svg.layers{1}
% svg.layers{3}.polys{4}

edit=1;
[h,svg]=plotSVG(svg,edit);

saveSVG('./building2.svg',svg);

 
 

