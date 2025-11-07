function newmap = bluewhitered(m, white_number)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.

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

bottom_blue = [0 0 0.5];
botmiddle_blue = [0 0.5 1];
middle_white = [1 1 1];
topmiddle_red = [1 0 0];
top_red = [0.5 0 0];

% Find limits
lims = get(gca, 'CLim');

bluemap_basic = [bottom_blue; botmiddle_blue; middle_white];
redmap_basic = [middle_white; topmiddle_red; top_red];

bluepart = white_number - lims(1);
redpart = lims(2) - white_number;

if bluepart > 0 && redpart > 0
    % It has both blue and red
    ratio_blue = bluepart / (bluepart + redpart);
    bluelen = round(m*ratio_blue);
    redlen = m - bluelen;

    bluemap = resize_cmap(bluemap_basic,bluelen);
    redmap = resize_cmap(redmap_basic,redlen);
    newmap = [bluemap; redmap];

elseif bluepart <= 0
    % Just red (positive)
    redmap_basic = [middle_white; topmiddle_red; top_red];
    newmap = resize_cmap(redmap_basic,m);

else
    % Just blue (negative)
    bluemap_basic = [bottom_blue; botmiddle_blue; middle_white];
    newmap = resize_cmap(bluemap_basic,m);

end
end

function newmap = resize_cmap(oldmap,newlen)

    oldlen = length(oldmap);
    oldsteps = linspace(0, 1, oldlen);
    newsteps = linspace(0, 1, newlen);
    newmap = zeros(newlen, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        temp = interp1(oldsteps, oldmap(:,i), newsteps);
        newmap(:,i) = min(max(temp', 0), 1);
    end
end
