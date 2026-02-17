function arrow(p1,p2)
% arrow from p1 to p2: see https://www.mathworks.com/matlabcentral/answers/160487-how-can-i-draw-a-line-with-arrow-head-between-2-data-points-in-a-plot#answer_365407
    xapf = @(x,pos,xl) pos(3)*(x-min(xl))/diff(xl)+pos(1);                % 'x' Annotation Position Function
    yapf = @(y,pos,yl) pos(4)*(y-min(yl))/diff(yl)+pos(2);                % 'y' Annotation Position Function
    xl = xlim;
    yl = ylim;
    pos = gca().Position;
    annotation('arrow', xapf([p1(1) p2(1)],pos,xl), yapf([p1(2) p2(2)],pos,yl))
end