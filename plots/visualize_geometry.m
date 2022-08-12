% close all
% p = parameters();
% b = var_bounds(p);
% x = [b.D_sft_nom, b.D_i_ratio_nom, b.D_or_ratio_nom, b.M_nom, b.N_WEC_nom, b.D_int_nom, b.w_n_nom];
% visualize_geometry(x,p)

function visualize_geometry(x,p,mini,color)

if nargin<3
    mini = false; % whether to adjust formatting for a mini plot inset
end
if nargin<4
    color = {'b','k','r'};
    compare = false;
else
    compare = true;
    color = {color,color,color};
end

[D_f,D_s_ratio,h_f_ratio,T_s_ratio,~,~,~,~] = deal(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));

if ~mini && ~compare
    figure
end
D_s = D_s_ratio * D_f;
h_f = h_f_ratio * D_f;
% Geometric similarity float submergence
T_f = p.T_f_over_h_f * h_f;
% Geometric similarity to maintain constant damping ratio
% D_s sets D_d, T_s, h_d
D_d = p.D_d_over_D_s * D_s;
T_s = p.T_s_over_D_s * D_s;
h_d = p.h_d_over_D_s * D_s;
% Another ratio defined by design variable
h_s = 1/T_s_ratio * T_s;

% waves
x = linspace(-30,30,100);
Hs = 1.5;
T = 7.5;
hold on
waves = plot(x,Hs*cos(x*2*pi/T),'c');

% WEC
center_rect([0 (h_s/2-T_s) D_s h_s],color{1})       % vert col
center_rect([0 -T_s D_d h_d],color{2})              % reaction plate
center_rect([0 (h_f/2-T_f) D_f h_f],color{3})       % float

% for legend
plot(NaN,NaN,color{1})

axis equal
if ~mini
    improvePlot
else
    set(gcf, 'Color', 'white');
end
grid on
ylim([-50 10])
xlim([-20 20])
set(waves,'HandleVisibility','off')

end

function center_rect(vec,color)

[x_mid,y_mid,w,h] = deal(vec(1),vec(2),vec(3),vec(4));
x_left = x_mid - w/2;
y_botm = y_mid - h/2;
pos = [x_left y_botm w h];

rectangle('Position',pos,'LineWidth',3,'EdgeColor',color)

end


