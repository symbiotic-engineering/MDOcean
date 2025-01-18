clc; clear; close all;
%% PARAMETERS OF SIMULATION

% Both:
% a1_num = .5;
% a2_num = 1;
% d1_num = .5;
% d2_num = .25;
% h_num = 1.05;
% m0_num = 1;

% MEEM:
% N=M=K=10

% BEM:
% resolution=0.03

%% MEEM
total_MEEM_time = 0.205;
compute_and_plot_percentage = 50.1; %percent
compute_eigen_hydro_coeff_percentage = 49.9; %percent
compute_and_plot = total_MEEM_time*compute_and_plot_percentage/100;
compute_eigen_hydro_coeff = total_MEEM_time*compute_eigen_hydro_coeff_percentage/100;

%% BEM

% 1ST RUN (INCLUDING COMPILING) AND 2ND RUN
total_BEM_time_1n2 = 17.45;
solve_1n2 = 11.8;
init_1n2 = 2.71;
mesh_vertical_cylinder_1n2 = 2.18;
immersed_part_1n2 = 542.65e-3;
add_1n2 = 297.52e-3;
face_centers_1n2 = 16.60e-3;

% 1ST RUN (INCLUDING COMPILING)
total_BEM_time_1 = 14.12;
solve_1 = 11.10;
init_1 = 1.45;
mesh_vertical_cylinder_1 = 1.13;
immersed_part_1 = 274.71e-3;
add_1 = 197.93e-3;
face_centers_1 = 13.81e-3;

% 2ND RUN ONLY
% total_BEM_time = total_BEM_time_1n2-total_BEM_time_1;
solve = solve_1n2-solve_1;
init = init_1n2-init_1;
mesh_vertical_cylinder = mesh_vertical_cylinder_1n2-mesh_vertical_cylinder_1;
immersed_part = immersed_part_1n2-immersed_part_1;
add = add_1n2-add_1;
face_centers = face_centers_1n2-face_centers_1;
total_BEM_time = (solve+init+...
    mesh_vertical_cylinder+immersed_part+add+face_centers);

%% NORMALIZE

% MEEM
compute_and_plot_norm = compute_and_plot/total_BEM_time;
compute_eigen_hydro_coeff_norm = compute_eigen_hydro_coeff/total_BEM_time;

% BEM
solve_n = solve/total_BEM_time;
init_n = init/total_BEM_time;
mesh_vertical_cylinder_n = mesh_vertical_cylinder/total_BEM_time;
immersed_part_n = immersed_part/total_BEM_time;
add_n = add/total_BEM_time;
face_centers_n = face_centers/total_BEM_time;
% other_n = other/total_BEM_time;

%% PLOT

legend_font_size = 12;
y_max = 1.01;

figure (1)
tiledlayout(2,1,"TileSpacing","none")

nexttile
x_MEEM = "MEEM";
y_MEEM = [compute_and_plot_norm compute_eigen_hydro_coeff_norm];
barh(x_MEEM,y_MEEM,'stacked')
xlim([0 y_max])
set(gca,'XGrid','on','YGrid','off')
ax = gca; 
ax.FontSize = 14;
hold on
plot([0.015, 0.03],[1, 1.75],'k-')
text(0.04,1.75,'Hydro Coefficients', 'FontSize',...
    legend_font_size,'Interpreter','latex')
plot([0.045, 0.06],[1, 0.25],'k-')
text(0.07,0.25,'Eigen Coefficients', 'FontSize',...
    legend_font_size,'Interpreter','latex')
ax.TickLabelInterpreter='latex';
x_axis_vec = 0:0.1:1;
set(gca, 'XTick', x_axis_vec, 'XTickLabel', string(x_axis_vec))
set(gca,'Xticklabel',[])


nexttile
x_BEM = "BEM";
y_BEM = [0 face_centers_n add_n immersed_part_n...
    solve_n mesh_vertical_cylinder_n init_n];
barh(x_BEM,y_BEM,'stacked')
xlim([0 y_max])
% set(gca,'Xticklabel',[])
set(gca,'XGrid','on','YGrid','off')
ax = gca; 
ax.FontSize = 14;
hold on
plot([0.015, 0.03],[1, 1.75],'k-')
text(0.04,1.75,'Add', 'FontSize', legend_font_size,'Interpreter','latex')
plot([0.07, 0.085],[1, 0.25],'k-')
text(0.095,0.25,'Immersed Part', 'FontSize',...
    legend_font_size,'Interpreter','latex')
plot([0.215, 0.215+0.015],[1, 1.75],'k-')
text(0.215+0.025,1.75,'Solve', 'FontSize',...
    legend_font_size,'Interpreter','latex')
plot([0.475, 0.475+0.015],[1, 0.25],'k-')
text(0.475+0.025,0.25,'Mesh Vertical Cylinder', 'FontSize',...
    legend_font_size,'Interpreter','latex')
plot([0.815, 0.815+0.015],[1, 1.75],'k-')
text(0.815+0.025,1.75,'Init', 'FontSize',...
    legend_font_size,'Interpreter','latex')
set(gca, 'XTick', x_axis_vec, 'XTickLabel', string(x_axis_vec))
xlabel('Normalized Run Time', 'FontSize', 14)
ax.TickLabelInterpreter='latex';
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

