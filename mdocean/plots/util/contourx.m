function [C,h] = contourx(X,Y,Z,varargin)
% wrapper around contourf that avoids "contour not rendered for constant
% zdata" warning by using imagesc when Z is constant

    Z_is_constant = numel(unique(Z(~isnan(Z)))) == 1;

    if Z_is_constant
        x = [min(X,[],'all') max(X,[],'all')];
        y = [min(Y,[],'all') max(Y,[],'all')];
        h = imagesc('XData',x,'YData',y,'CData',Z,'AlphaData',~isnan(Z));
        C = [];
    else
        [C,h] = contourf(X,Y,Z,varargin{:});
    end
end
