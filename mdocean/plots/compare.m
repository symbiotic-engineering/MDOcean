function [DV_table,out_table] = compare(filename_uuid)

if nargin==0
    filename_uuid = '';
end

p = parameters();
b = var_bounds();
b.filename_uuid = filename_uuid;
x0_input = b.X_start_struct;

[Xs_opt,~,~,~,~,~,~,vals_opt] = gradient_optim(x0_input,p,b);
X_minLCOE = Xs_opt(:,1);
X_minCapex = Xs_opt(:,2);
val_minLCOE = vals_opt(1);
val_minCapex = vals_opt(2);

[X_maxPower,val_maxPower] = max_avg_power(p,b);

p_bal = p;
p_bal.avg_power_min = b.power_balanced;
[X_balanced,~,~,~,~,~,~,val_balanced] = gradient_optim(x0_input,p_bal,b,2)

X_nom	= [b.X_noms' 1];
[~,~,~,~,val_nom] = simulation(X_nom,p);

%%
X    = [X_nom;   X_minLCOE';  X_minCapex';  X_maxPower';  X_balanced'];
vals = [val_nom, val_minLCOE, val_minCapex, val_maxPower, val_balanced];

titles = {'Nominal','Min LCOE','Min Capex','Max Power','Balanced'};
color = {'k','b','r','g','m'};
num_designs = length(titles);

%% geometry comparison
figure
% sort by float diameter
idx_D_f = strcmp(b.var_names,'D_f');
D_f = X(:,idx_D_f);
[~,biggest_to_smallest] = sort(D_f,'descend');

% find x range
idx_D_s = strcmp(b.var_names,'D_s');
D_d = p.D_d_over_D_s * X(:,idx_D_s);
biggest_diam = max([D_f;D_d]);
x_max = biggest_diam/2 * 1.1;

% find y range
T_s = p.T_s_over_D_s * X(:,idx_D_s);
y_min = -max(T_s) * 1.1;

% waves
x = linspace(-x_max,x_max,100);
Hs = 3;
T = 7.5;
hold on
waves = plot(x,Hs/2*cos(x*2*pi/T),'c','Linewidth',3);
set(waves,'HandleVisibility','off')

for i=1:num_designs
    x = X(biggest_to_smallest(i),:);
    hold on
    visualize_geometry(x,p,false,color{biggest_to_smallest(i)})
end
xlim([-x_max,x_max])
y_max = max([findobj(gca().Children,'Type','Line').YData]) * 1.1;
ylim([y_min y_max])
legend(titles(biggest_to_smallest),'Position',[0.5941    0.3508    0.2663    0.2167])

%% power probability comparison
f = figure;
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
ylabel(t,'Probability (-)','FontWeight','bold')
for i=1:num_designs
    x = X(i,:);
    hold on
    t = power_PDF(x,p,color{i},t,[false true]);
end
legend(titles{:},'location','best')
delete(nexttile(1)) % delete blank PDF so just CDF remains
f.Position(3:4) = [792 484];

%% power matrix comparison

figure
for i=1:num_designs
    x = X(i,:);
    [~,~,P_matrix] = simulation(x, p);
    P_matrix = P_matrix / 1e3; % convert W to kW
    [T,Hs] = meshgrid(p.T,p.Hs);

    subplot(2,3,i);
    hold on
    contourf(T,Hs,P_matrix .* p.JPD/100)
    improvePlot
    grid on
    title(titles{i})
    caxis([0 15])
    ax = gca;
    second_col = i==2 || i==5;
    third_col = i==3;
    if second_col
        xoffset = -ax.Position(3)*.25;
    elseif third_col
        xoffset = -ax.Position(3)*.5;
    else
        xoffset = 0;
    end
    second_row = i>3;
    if second_row
        yoffset = ax.Position(4)*.2;
    else
        yoffset = 0;
    end
    set(ax,'Position',[ax.Position(1)+xoffset ax.Position(2)+yoffset ...
        ax.Position(3)*.8 ax.Position(4)*.8])
    set(ax.Title,'FontWeight','normal')
    set(ax,'Layer','top','GridColor',[.9 .9 .9])
end
cb = colorbar('Position',[.8 .175 .08 .685]);
%cb.Label.String = 'Weighted Power (kW)';
sgtitle('Weighted Power (kW)','FontWeight','bold','FontSize',20)
text(-10,-2,'Wave Period T (s)','FontWeight','bold','FontSize',16)
text(-30,5,'Wave Height Hs (m)','FontWeight','bold','FontSize',16,'Rotation',90)

%% hydro coeff comparison plot
hydro_compare(vals,color)

%% design variable table
DV_table = array2table(X.', ...
    'VariableNames',titles, 'RowNames', b.var_names_pretty);

%% output table
temp_table = struct2table(vals,'RowNames',titles);
scalar_vals = varfun(@isnumeric, temp_table, 'OutputFormat', 'uniform');
hydro_names = {'over_rho','phase'};
hydro_coeff_rows = contains(temp_table.Properties.VariableNames,hydro_names);
out_table = rows2vars(temp_table(:,scalar_vals & ~hydro_coeff_rows));
out_table.Properties.RowNames = out_table.OriginalVariableNames;
out_table = removevars(out_table,'OriginalVariableNames');

end

function hydro_compare(vals,colors)
figure
subplot(1,2,1)
hold on
dummy_style = {['k' '-'],['b' '-'],['r' '-'],['g' '-'],['m' '-'],...
    ['k' '-*'],['k' '-.']};
for i = 1:length(dummy_style)
    plot(NaN,NaN,dummy_style{i})
end
colororder({'k','k'})
for i=1:length(vals)
    val = vals(i);
    col = colors{i};
    w=unique(val.w(~isnan(val.w)),'stable');
    gamma_f_over_rho_g = unique(val.gamma_f_over_rho_g(~isnan(val.gamma_f_over_rho_g)),'stable');
    gamma_phase_f = unique(val.gamma_phase_f(~isnan(val.gamma_phase_f)),'stable');
    yyaxis left
    plot(w(:,1), gamma_f_over_rho_g(:,1),[col '-*'])
    ylabel('\gamma_{f}/(\rhog) (m^{2})')
    yyaxis right
    plot(w(:,1), gamma_phase_f(:,1),     [col '-.'])
    ylabel('\angle\gamma_{f} (rad)')
end
title('Hydrodynamic Coefficients')
xlabel('Wave Frequency (\omega)')
xlim([0.3,1.45])
leg = legend({'Nominal','Min LCOE','Min CAPEX','Max Power','Balanced', ...
    '\gamma_{f}/(\rhog)','\angle\gamma_{f}'});
hold off
improvePlot
leg.Location='northeast';


subplot(1,2,2)
hold on
dummy_style = {['k' '-'],['b' '-'],['r' '-'],['g' '-'],['m' '-'],...
    ['k' '--'],['k' '-']};
for i = 1:length(dummy_style)
    plot(NaN,NaN,dummy_style{i})
end
for i=1:length(vals)
    val = vals(i);
    col = colors{i};
    w=unique(val.w(~isnan(val.w)),'stable');
    A_f_over_rho = unique(val.A_f_over_rho(~isnan(val.A_f_over_rho)),'stable');
    B_f_over_rho_w = unique(val.B_f_over_rho_w(~isnan(val.B_f_over_rho_w)),'stable');
    plot(w(:,1), A_f_over_rho(:,1),[col '--'],...
        w(:,1), B_f_over_rho_w(:,1),     [col '-'])
end
title('Hydrodynamic Coefficients')
xlabel('Wave Frequency (\omega)')
xlim([0.3,1.45])
%ylim([-10, 26900])
leg2=legend({'Nominal','Min LCOE','Min CAPEX','Max Power','Balanced', ...
    'A_{f}/\rho','B_{f}/(\rho\omega)'});
ylabel('A_{f}/\rho (m^{3}), B_{f}/(\rho\omega) (m^{3}/(rad^{2}s^{2}))')
hold off
improvePlot
leg2.Location='northeast';
h=findobj(gca().Children,"Type","line");
end