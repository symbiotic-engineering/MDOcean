function [DV_table,out_table] = compare(filename_uuid)

if nargin==0
    filename_uuid = '';
end

p = parameters();
b = var_bounds();
b.filename_uuid = filename_uuid;
x0_input = b.X_start_struct;

[Xs_opt,~,~,~,vals_opt] = gradient_optim(x0_input,p,b);
X_minLCOE = Xs_opt(:,1);
X_minCapex = Xs_opt(:,2);
val_minLCOE = vals_opt(1);
val_minCapex = vals_opt(2);

[X_maxPower,val_maxPower] = max_avg_power(p,b);

p.LCOE_max = 0.1;
[X_balanced,~,~,~,val_balanced] = gradient_optim(x0_input,p,b,2);

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
biggest_to_smallest = [5, 2, 4, 1, 3];
for i=1:num_designs
    x = X(biggest_to_smallest(i),:);
    hold on
    visualize_geometry(x,p,false,color{i})
end
legend({'Balanced','Min LCOE','Max Power','Nominal','Min CAPEX'},'Location','east')

%% power probability comparison
figure
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
ylabel(t,'Probability (-)','FontWeight','bold')
for i=1:num_designs
    x = X(i,:);
    hold on
    t = power_PDF(x,p,color{i},t);
end
legend(titles{:},'location','best')

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

%% design variable table
DV_table = array2table(X.', ...
        'VariableNames',titles, 'RowNames', b.var_names_pretty);

%% output table
temp_table = struct2table(vals,'RowNames',titles);
scalar_vals = varfun(@isnumeric, temp_table, 'OutputFormat', 'uniform');
out_table = rows2vars(temp_table(:,scalar_vals));
out_table.Properties.RowNames = out_table.OriginalVariableNames;
out_table = removevars(out_table,'OriginalVariableNames');

end