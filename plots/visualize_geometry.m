% close all
% p = parameters();
% b = var_bounds(p);
% x = [b.D_sft_nom, b.D_i_ratio_nom, b.D_or_nom, b.M_nom, b.N_WEC_nom, b.D_int_nom, b.w_n_nom];
% visualize_geometry(x,p)

function visualize_geometry(x,p,mini)

if nargin<3
    mini = false; % whether to adjust formatting for a mini plot inset
end

[D_sft,D_i_ratio,D_or,~,N_WEC,D_int,w_n] = deal(x(1),x(2),x(3),x(4),x(5),x(6),x(7));

if ~mini
    figure
end
t_f = D_sft / 4;
D_i = D_i_ratio * D_sft;
h = D_i * 7;

% waves
x = linspace(-30,30,100);
Hs = 1.5;
T = 7.5;
hold on
plot(x,Hs*cos(x*2*pi/T),'c')

% WEC
center_rect([0 (t_f-h)/2+5 D_i h],'b')      % vert col
center_rect([0 -h+5+t_f/2 D_or p.t_r],'k')  % reaction plate
center_rect([0 0 D_sft t_f],'r')            % float

axis equal
if ~mini
    improvePlot
else
    set(gcf, 'Color', 'white');
end
grid on
ylim([-60 10])
xlim([-30 30])

end

function center_rect(vec,color)

[x_mid,y_mid,w,h] = deal(vec(1),vec(2),vec(3),vec(4));
x_left = x_mid - w/2;
y_botm = y_mid - h/2;
pos = [x_left y_botm w h];

rectangle('Position',pos,'LineWidth',3,'EdgeColor',color)

end


