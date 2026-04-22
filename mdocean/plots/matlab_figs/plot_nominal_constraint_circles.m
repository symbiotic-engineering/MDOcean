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
if nargin>1
    existing_handles = existing_axes.Children;
    ax = axes;
    new_handles = copyobj(existing_handles,ax);
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
axis equal;
pbaspect([1 1 1]);

% Draw circles and shade feasible interiors
theta_vec = linspace(0, 2*pi, 360);
colors = lines(N);

for k = 1:N
    cx = centers(k,1);  cy = centers(k,2);  r = radii(k);
    px = cx + r * cos(theta_vec);
    py = cy + r * sin(theta_vec);
    fill(px, py, colors(k,:), 'FaceAlpha', 0.15, 'EdgeColor', colors(k,:), ...
         'LineWidth', 1.5, 'DisplayName', labels{k});
end

% Mark origin
plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Origin (unconstrained opt)');

% Mark constrained optimal
if isfinite(Gamma_opt)
    plot(real(Gamma_opt), imag(Gamma_opt), 'gp', 'MarkerFaceColor', 'g', ...
         'MarkerSize', 16, 'DisplayName', 'Constrained optimal \Gamma');
end

% Draw arrows around each circle pointing inward (toward center)
n_arrows = 6;
arrow_len = min(radii) * 0.25;
for k = 1:N
    cx = centers(k,1);  cy = centers(k,2);  r = radii(k);
    for m = 1:n_arrows
        ang = (m-1) * 2*pi / n_arrows;
        px = cx + r * cos(ang);  py = cy + r * sin(ang);
        dx = cx - px;  dy = cy - py;
        len = sqrt(dx^2 + dy^2);
        p1 = [px, py];
        p2 = [px + arrow_len * dx/len, py + arrow_len * dy/len];
        arrow(p1, p2);
    end
end

% Add text labels near circle boundaries at top of each circle
for k = 1:N
    cx = centers(k,1);  cy = centers(k,2);  r = radii(k);
    tx = cx;  ty = cy + r + 0.05 * r;
    text(tx, ty, labels{k}, 'HorizontalAlignment', 'center', ...
         'Color', colors(k,:), 'FontSize', 9, 'FontWeight', 'bold');
end

% Add period info in title
T_i = 2*pi / w_i;
title(sprintf('QCQP constraint circles (T = %.1f s)', T_i));
xlabel('Re(\Gamma)');
ylabel('Im(\Gamma)');
legend('Location', 'eastoutside');
improvePlot;
end
