function newmap = bluewhitered(m, white_number, flip, full_color)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
% white_number (default zero) is the number to correspond to white.
% (Sometimes this is set to 1 if plotting a ratio).
% flip (default false) is a boolean for whether to use red-white-blue instead
% of blue-white-red.
% full_color (default false) is a boolean for whether the entire range of
% blue and red shades should be shown even if it results in uneven color
% gradient. for example, numerical limits of [-1,100] will have colors of
% [nearly-white-pale-blue, full-red] when full_color is false and
% [full-blue, full-red] when full_color is true.

% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gca,'colormap'),1);
end
if nargin < 2
    white_number = 0;
end
if nargin < 3
    flip = false;
end
if nargin < 4
    full_color = false;
end

% define rgb colors
blue       = [0   0   0.5];
light_blue = [0   0.5 1];
white      = [1   1   1];
red        = [1   0   0];
dark_red   = [0.5 0   0];

% Find limits
lims = get(gca, 'CLim');

if ~flip % blue white red
    botmap_basic = [blue; light_blue; white];
    topmap_basic = [white; red; dark_red];
else % red white blue
    botmap_basic = [dark_red; red; white];
    topmap_basic = [white; light_blue; blue];
end

botpart = white_number - lims(1);
toppart = lims(2) - white_number;

if botpart > 0 && toppart > 0
    % It has both bot (blue) and top (red)
    ratio_bot = botpart / (botpart + toppart);
    botlen = round(m*ratio_bot);
    toplen = m - botlen;
    
    if full_color || botlen == toplen
        frac_bot = 0;
        frac_top = 1;
    elseif botlen > toplen
        frac_bot = 0;
        frac_top = toplen/botlen;
    else
        frac_bot = 1-botlen/toplen;
        frac_top = 1;
    end

    botmap = resize_cmap(botmap_basic,botlen,frac_bot,1);
    topmap = resize_cmap(topmap_basic,toplen,0,frac_top);
    newmap = [botmap; topmap];
    
elseif botpart <= 0
    % Just top (red, positive)
    frac = -botpart/toppart;
    newmap = resize_cmap(topmap_basic,m,frac,1);
    
else
    % Just bottom (blue, negative)
    frac = 1+toppart/botpart;
    newmap = resize_cmap(botmap_basic,m,0,frac);
    
end
end

function newmap = resize_cmap(oldmap,newlen,newstep_min,newstep_max)

    if nargin==2
        newstep_min = 0;
        newstep_max = 1;
    end

    oldlen = length(oldmap);
    oldsteps = linspace(0, 1, oldlen);
    newsteps = linspace(newstep_min, newstep_max, newlen);
    newmap = zeros(newlen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        temp = interp1(oldsteps, oldmap(:,i), newsteps);
        newmap(:,i) = min(max(temp', 0), 1);
    end
end
