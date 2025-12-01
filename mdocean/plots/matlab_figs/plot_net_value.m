clear;close all

p = parameters();
in = p;
b = var_bounds();
X = [b.X_noms;1];
[~, P_matrix_elec, ~, val] = simulation(X, p);

%calculate revenue and avoided carbon
%% power/energy
% power
P_weighted = P_matrix_elec .* p.JPD / 100 * in.eff_array; % taken from dynamics.m and LCOE_from_capex_design_power.m
P_avg = sum(P_weighted,'all','omitnan');

% annual energy
hr_per_yr = 8766;
AEP_matrix = P_weighted * in.N_WEC * hr_per_yr / 1000; % W to kWh per year, all wecs
AEP = P_avg * in.N_WEC * hr_per_yr / 1000; % kWh per year

%% force for scaling
%scale the cost by the amount of steel needed at a given sea state
%capex equals the sum of the entries of capex_matrix
factor = 1./sum(val.F_heave_mat,'all','omitnan');
force_matrix = val.F_heave_mat.*factor;

%% econ

levelized_cost = in.FCR * val.capex + val.opex; % $/year
levelized_cost_matrix_unscaled = in.FCR * val.capex * force_matrix + val.opex; % $/year
cost_matrix_scale_factor = levelized_cost / mysum(levelized_cost_matrix_unscaled);
levelized_cost_matrix = levelized_cost_matrix_unscaled * cost_matrix_scale_factor;

levelized_value_matrix = AEP_matrix .* p.marginal_price; % kWh per year * $/kWh = $/year

[LCOE_matrix,LCOE,...
 LVOE_matrix,LVOE,...
 NVOE_matrix,NVOE,...
 net_value_matrix,net_value,...
 VCR_matrix, VCR] = get_LXOE_NVOE_VCR(levelized_cost_matrix, ...
                                      levelized_value_matrix, ...
                                      AEP_matrix);

% lcoe check
assert(ismembertol(val.LCOE,LCOE));

% final var renaming
econ_cost = LCOE_matrix;
econ_value = LVOE_matrix; 

%% enviro

levelized_eco_value_matrix = AEP_matrix .* p.marginal_carbon * p.eco_cost_carbon; % kWh/year * kgCO2/kWh * $/kgCO2 = $/year
levelized_eco_cost = val.env_cost_per_wec * p.N_WEC; % fixme this isn't levelized
levelized_eco_cost_matrix = levelized_eco_cost * force_matrix;

[LCOE_matrix_eco,LCOE_eco,...
 LVOE_matrix_eco,LVOE_eco,...
 NVOE_matrix_eco,NVOE_eco,...
 net_value_matrix_eco,net_value_eco,...
 VCR_matrix_eco, VCR_eco] = get_LXOE_NVOE_VCR(levelized_eco_cost_matrix, ...
                                              levelized_eco_value_matrix, ...
                                              AEP_matrix);

%% wave conditions
[T_mdocean,Hs_mdocean] = meshgrid(p.T,p.Hs);

%% figures
%plot net economic value multiplication
f1 = figure;
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
contourf(Hs_mdocean, T_mdocean, levelized_cost_matrix)
title("Levelized cost (scaled by force)")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'=','FontSize',50)
axis off

ax1 = nexttile;
contourf(ax1, Hs_mdocean, T_mdocean, VCR_matrix)
colormap(ax1, bluewhitered(256,1,true))
title("Economic Value to Cost Ratio")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar(ax1)
improvePlot

%plot net environmental value multiplication
f2 = figure;
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
contourf(Hs_mdocean, T_mdocean, levelized_eco_cost_matrix)
title("Environmental Cost")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar

nexttile
text(0.5,0.5,'=','FontSize',50)
axis off

ax2 = nexttile;
contourf(ax2, Hs_mdocean, T_mdocean, VCR_matrix_eco)
colormap(ax2, bluewhitered(256,1,true))
title("Environmental Value to Cost Ratio")
xlabel("Wave Height")
ylabel("Wave Period")
colorbar(ax2)
improvePlot

%plot net value
f3 = figure;
tiledlayout(1,2)
ax1 = nexttile;
contourf(ax1, Hs_mdocean, T_mdocean, net_value_matrix)
colormap(ax1, bluewhitered(256,0,true))
title("Net Economic Value")
colorbar(ax1)

ax2 = nexttile;
contourf(ax2, Hs_mdocean, T_mdocean, net_value_matrix_eco)
colormap(ax2, bluewhitered(256,0,true))
title("Net Environmental Value")
colorbar(ax2)
improvePlot


function [LCOE_matrix,LCOE,...
          LVOE_matrix,LVOE,...
          NVOE_matrix,NVOE,...
          net_value_matrix,net_value,...
          VCR_matrix, VCR] = get_LXOE_NVOE_VCR(levelized_cost_matrix, ...
                                               levelized_value_matrix, ...
                                               AEP_matrix)
    
    % lcoe = levelized cost / annual energy production
    [LCOE_matrix, LCOE] = LXOE(levelized_cost_matrix,  AEP_matrix);

    % lvoe = levelized value / annual energy production
    [LVOE_matrix, LVOE] = LXOE(levelized_value_matrix, AEP_matrix);

    % net value, NVOE, and value cost ratio matrices
    net_value_matrix = levelized_value_matrix - levelized_cost_matrix;
    NVOE_matrix      = LVOE_matrix            -  LCOE_matrix;
    VCR_matrix       = LVOE_matrix           ./  LCOE_matrix;

    % net value, NVOE, and value cost ratio scalars
    net_value = mysum(net_value_matrix);
    NVOE = LVOE - LCOE;
    VCR  = LVOE / LCOE;
    
end

function [LXOE_matrix,LXOE] = LXOE(levelized_X_matrix,AEP_matrix)
    % lxoe = levelized X (where X = cost or value) / annual energy production
    LXOE_matrix = levelized_X_matrix ./ AEP_matrix; % from LCOE_from_capex_design_power.m: $/year / kWh/year = $/kWh
    levelized_X = mysum(levelized_X_matrix);
    AEP = mysum(AEP_matrix);
    LXOE = levelized_X ./ AEP;
end

function out = mysum(in)
    out = sum(in,'all','omitnan');
end