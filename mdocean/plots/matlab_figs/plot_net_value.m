p=parameters();
in=p;
b=var_bounds();
X=[b.X_noms;1];
[~, P_matrix_elec, ~, val] = simulation(X, p);

%calculate revenue and avoided carbon
econ_value = P_matrix_elec .* p.marginal_price .* p.JPD;
env_value = P_matrix_elec .* p.marginal_carbon .* p.JPD;
net_econ_value = econ_value ./ val.capex;
net_env_value = env_value ./ val.env_cost_per_wec;
[T_mdocean,Hs_mdocean] = meshgrid(p.T,p.Hs);

%scale the cost by the amount of steel needed at a given sea state
%capex equals the sum of the entries of capex_matrix
factor = 1./sum(val.F_heave_mat,'all','omitnan');
capex_matrix = val.capex*val.F_heave_mat.*factor;
env_cost_matrix = val.env_cost_per_wec*val.F_heave_mat.*factor;

%plot net economic value multiplication
figure
tiledlayout(2,5,'TileSpacing','compact','Padding','compact')
nexttile
contourf(Hs_mdocean, T_mdocean, P_matrix_elec)
title("Power")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'x','FontSize',50)
axis off

nexttile
contourf(Hs_mdocean, T_mdocean, in.marginal_price)
title("Marginal Price")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'x','FontSize',50)
axis off

nexttile
contourf(Hs_mdocean, T_mdocean, p.JPD)
title("JPD")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,char(247),'FontSize',50)
axis off

nexttile
contourf(Hs_mdocean, T_mdocean, capex_matrix)
title("Capex")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'=','FontSize',50)
axis off

ax1 = nexttile;
contourf(ax1, Hs_mdocean, T_mdocean, net_econ_value)
colormap(ax1, bluewhitered)
title("Economic Value to Cost Ratio")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar(ax1)
improvePlot

%plot net environmental value multiplication
figure
tiledlayout(2,5,'TileSpacing','compact','Padding','compact')
nexttile
contourf(Hs_mdocean, T_mdocean, P_matrix_elec)
title("Power")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'x','FontSize',50)
axis off

nexttile
contourf(Hs_mdocean, T_mdocean, in.marginal_carbon)
title("Marginal Carbon")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'x','FontSize',50)
axis off

nexttile
contourf(Hs_mdocean, T_mdocean, p.JPD)
title("JPD")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,char(247),'FontSize',50)
axis off

nexttile
contourf(Hs_mdocean, T_mdocean, env_cost_matrix)
title("Environmental Cost")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'=','FontSize',50)
axis off

ax2 = nexttile;
contourf(ax2, Hs_mdocean, T_mdocean, net_env_value)
colormap(ax2, bluewhitered)
title("Environmental Value to Cost Ratio")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar(ax2)
improvePlot

%plot net value
figure
tiledlayout(1,2)
ax1 = nexttile;
contourf(ax1, Hs_mdocean, T_mdocean, net_econ_value)
colormap(ax1, bluewhitered)
title("Net Economic Value")
colorbar(ax1)

ax2 = nexttile;
contourf(ax2, Hs_mdocean, T_mdocean, net_env_value)
colormap(ax2, bluewhitered)
title("Net Environmental Value")
colorbar(ax2)
improvePlot