% close all
% p = parameters();
% b = var_bounds();
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

[D_s,D_f,T_f_2,h_s,h_fs_clear,~,~,~] = deal(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));

if ~mini && ~compare
    figure
end

% Geometric similarity float submergence
T_f_1 = p.T_f_1_over_T_f_2 * T_f_2;
D_f_b = p.D_f_b_over_D_f * D_f;
h_f = T_f_2 / p.T_f_2_over_h_f;
% Geometric similarity to maintain constant damping ratio
% D_s sets D_d, T_s, h_d
D_d = p.D_d_over_D_s * D_s;
T_s = p.T_s_over_D_s * D_s;
h_d = p.h_d_over_D_s * D_s;

% waves
x = linspace(-30,30,100);
Hs = 3;
T = 7.5;
hold on
waves = plot(x,Hs/2*cos(x*2*pi/T),'c');

% WEC
float_rect_bottom = -T_f_1;
float_rect_height = h_f - (T_f_2-T_f_1);
float_rect_middle = float_rect_bottom + float_rect_height/2;

center_rect([0 (h_s/2-T_s) D_s h_s],color{1})       % vert col
center_rect([0 -T_s D_d h_d],color{2})              % reaction plate
center_rect([0 float_rect_middle D_f float_rect_height],color{3})       % float straight part

trapezoid(D_f_b,D_f,-T_f_2,-T_f_1,color{3}) % float slanted part

% float tubes
m_tubes = (h_fs_clear + (h_s-T_s) - (h_f - T_f_2)) / ( (D_f-D_s)/2 );
x_tubes = [-D_f/2, 0, D_f/2];
y_tube_1 = h_f - T_f_2;
y_tubes = y_tube_1 + [0, m_tubes * D_f/2, 0];
h = plot(x_tubes, y_tubes, color{3},'LineWidth',3);
h.Annotation.LegendInformation.IconDisplayStyle= 'off'; % no legend

% for legend
plot(NaN,NaN,color{1})

axis equal
if ~mini
    improvePlot
else
    set(gcf, 'Color', 'white');
end
grid on
ylim([-40 20])
xlim([-27 27])
set(waves,'HandleVisibility','off')

end

function center_rect(vec,color)

    [x_mid,y_mid,w,h] = deal(vec(1),vec(2),vec(3),vec(4));
    x_left = x_mid - w/2;
    y_botm = y_mid - h/2;
    pos = [x_left y_botm w h];
    
    rectangle('Position',pos,'LineWidth',3,'EdgeColor',color)

end

function trapezoid(base_1,base_2,y_1,y_2,color)
    x = [-base_1/2, base_1/2, base_2/2, -base_2/2, -base_1/2];
    y = [y_1 y_1 y_2 y_2 y_1];
    h = plot(x,y,'LineWidth',3,'Color',color);
    h.Annotation.LegendInformation.IconDisplayStyle= 'off'; % no legend
end
