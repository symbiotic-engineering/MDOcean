svg=loadSVG('building.svg');
svg
% svg.images{1}
% svg.layers{1}
% svg.layers{3}.polys{4}
pwd
edit=1
[h,svg]=plotSVG(svg,edit);pwd

saveSVG('./building2.svg',svg);

 
 

