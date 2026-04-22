function fig = plot_nominal_constraint_circles(val, existing_axes)
% PLOT_NOMINAL_CONSTRAINT_CIRCLES  Visualize QCQP constraint circles in Gamma
%   plane for the nominal design at the most constrained sea state.
%   Shows the force, amplitude, and positive-power constraint circles with
%   the feasible region shaded, directional arrows, and constraint labels.
%
%   Inputs:
%     val  - simulation output struct (must contain val.qcqp_debug)
%
%   Output:
%     fig  - figure handle

fig = figure;
pax = polaraxes;

if nargin>1
    existing_handles = existing_axes.Children;

    new_handles = copyobj(existing_handles,pax);
    colorbar
    hold on
end

if ~isfield(val, 'qcqp_debug') || isempty(val.qcqp_debug.centers)
    title('No constrained sea state found');
    return
end

d = val.qcqp_debug;
centers = d.centers;        % N×2 array of [Re Im] circle centers
radii   = d.radii;          % N×1 vector of radii
labels  = d.labels;         % N×1 cell of constraint names
Gamma_opt = d.Gamma_opt;    % complex optimal Gamma
w_i     = d.w;              % angular frequency of this sea state

N = size(centers, 1);

hold on;
%axis equal;
%pbaspect([1 1 1]);

% Draw circles and shade feasible interiors
colors = lines(N);
transparency = 0.15;

cax = fillcircles(centers,radii,pax,colors,transparency,labels);

hold on
% Mark origin
plot(pax, 0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Origin (unconstrained opt)');

% Mark constrained optimal
if isfinite(Gamma_opt)
    plot(pax,angle(Gamma_opt), abs(Gamma_opt), 'gp', 'MarkerFaceColor', 'g', ...
         'MarkerSize', 16, 'DisplayName', 'Constrained optimal \Gamma');
end

% Draw arrows around each circle pointing inward (toward center)
% n_arrows = 6;
% arrow_len = min(radii) * 0.25;
% for k = 1:N
%     cx = centers(k,1);  cy = centers(k,2);  r = radii(k);
%     for m = 1:n_arrows
%         ang = (m-1) * 2*pi / n_arrows;
%         px = cx + r * cos(ang);  py = cy + r * sin(ang);
%         dx = cx - px;  dy = cy - py;
%         len = sqrt(dx^2 + dy^2);
%         p1 = [px, py];
%         p2 = [px + arrow_len * dx/len, py + arrow_len * dy/len];
%         arrow(p1, p2);
%     end
% end

% Add text labels near circle boundaries at top of each circle
for k = 1:N
    cx = centers(k,1);  cy = centers(k,2);  r = radii(k);
    tx = cx;  ty = cy - r - 0.05 * r;
    text(cax, tx, ty, labels{k}, 'HorizontalAlignment', 'center', ...
         'Color', colors(k,:), 'FontSize', 9, 'FontWeight', 'bold');
end

% Add period info in title
T_i = 2*pi / w_i;
title(pax, sprintf('QCQP constraint circles (T = %.1f s)', T_i));
xlabel(cax, '$\Re(\Gamma)$','Interpreter','latex');
ylabel(cax, '$\Im(\Gamma)$','Interpreter','latex');
cax.XAxis.Label.Color=[0 0 0];
cax.XAxis.Label.Visible='on';
cax.YAxis.Label.Color=[0 0 0];
cax.YAxis.Label.Visible='on';
improvePlot;
end

function cartesianAx = fillcircles(centers,radii,polarAx,colors,transparency,labels)
% adapted from https://www.mathworks.com/matlabcentral/answers/1696275-how-to-make-a-fill-for-standard-deviation-around-a-polar-plot

    % Overlay a new Cartesian axes on the given polar axes
    cartesianAx = axes('Position', polarAx.Position, 'Color', 'none');
%     cartesianAx.XAxis.Visible = 'off';
%     cartesianAx.YAxis.Visible = 'off';
    hold(cartesianAx,'on')

    radiusLimit = max(get(polarAx, 'RLim'));
    theta_vec = linspace(0, 2*pi, 360);

    N = size(centers, 1);
    for k = 1:N
        cx = centers(k,1);  cy = centers(k,2);  r = radii(k);
        px = cx + r * cos(theta_vec);
        py = cy + r * sin(theta_vec);
        [ptheta,prho] = cart2pol(px,py);
        [px_cropped, py_cropped] = pol2cart(ptheta, min(prho,radiusLimit));
    
        fill(cartesianAx, px_cropped, py_cropped, colors(k,:), ...
         'FaceAlpha', transparency, 'EdgeAlpha', 0, ...
         'EdgeColor', colors(k,:), ...
         'LineWidth', 1.5, 'DisplayName', labels{k});
    end

    % Adjust limits of Cartesian axes to match polar plot and hide them

    legend(cartesianAx,'Location','best')

    xlim(cartesianAx, [-radiusLimit, radiusLimit]);
    ylim(cartesianAx, [-radiusLimit, radiusLimit]);
    axis(cartesianAx, 'square', 'off');
end


function cartesianAx = fillpolar(polarAx, cartesianAx, thetaLower, thetaUpper, radiusLower, radiusUpper,...
                    fillColor, transparency, label)

    % Transform polar coordinates to Cartesian for filling
    [xLower, yLower] = pol2cart(thetaLower, radiusLower);
    [xUpper, yUpper] = pol2cart(fliplr(thetaUpper), fliplr(radiusUpper));
    
    % Draw filled area between the curves using Cartesian coordinates
    
    
end
