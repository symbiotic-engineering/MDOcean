%   License FreeBSD:
%
%      Copyright (c) 2016  Martin de La Gorce
%      All rights reserved.
%
%      Redistribution and use in source and binary forms, with or without
%      modification, are permitted provided that the following conditions are met:
%
%      1. Redistributions of source code must retain the above copyright notice, this
%         list of conditions and the following disclaimer.
%      2. Redistributions in binary form must reproduce the above copyright notice,
%         this list of conditions and the following disclaimer in the documentation
%         and/or other materials provided with the distribution.
%
%      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%      ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%      ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%      (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%      LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%      ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%      SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%      The views and conclusions contained in the software and documentation are those
%      of the authors and should not be interpreted as representing official policies,
%      either expressed or implied, of the FreeBSD Project.


function [handles,svg]=plotSVG(svg,interactive)
if nargin<2
    interactive=false;
end

f=gcf();
set(f,'Renderer','OpenGL');
handles={};
for idImage=1: length(svg.images)
    image=svg.images{idImage};
    imagefile=image.file;
    handles{end+1}=imshow(imread(imagefile),'Xdata',[image.x,image.x+image.width],'Ydata',[image.y,image.y+image.height]);
    axis('ij');
    hold on;
end


for idLayer =1: length(svg.layers)
    layer=svg.layers{idLayer};
    polygons=layer.polys;
    colors=layer.colors;
    hold on;
    
    for k =1 : length(polygons)
        poly=polygons{k};
        if (~isempty(poly))
            h=fill(poly(1,:),poly(2,:),colors(1:3,k)');
            handles{end+1}=h;
            if size(colors,1)==4
                set(h,'facealpha',colors(4,k));
            end
            
            set(h,'UserData',[idLayer,k]);
        end
    end
    
    
end


selected=[];
idvertex=0;
if interactive
    set(gcf,'WindowButtonDownFcn',@down);
    set(gcf,'WindowButtonUpFcn',@up);
    set(gcf,'WindowButtonMotionFcn',@motion);
    drawnow();
    waitfor(f);
    fprintf('done')
    
end

    function down(~,~)
        if isempty(selected)
        
            p=get(gca,'currentpoint');
            p=p(1,1:2);
            selected=gco();
            if isa(selected,'image')
                selected=[];
            elseif not(isempty(selected))
               
                v=get(selected,'Vertices');
                d=sum((v-repmat(p,size(v,1),1)).^2,2);
                [ ~,idvertex]=min(d);
                
            end
        end
    end

    function up(~,~)
        selected=[];
      
    end

    function motion(~,~)
        if not(isempty(selected))
            p=get(gca,'currentpoint');
            p=p(1,1:2);
            selected=gco();
            v=get(selected,'Vertices');
            v(idvertex,:)=p;
            set(selected,'Vertices',v);
            d=get(selected,'UserData');
            svg.layers{d(1)}.polys{d(2)}=v';
        end
    end
end




