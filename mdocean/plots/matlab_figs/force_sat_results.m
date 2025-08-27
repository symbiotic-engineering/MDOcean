function [fig] = force_sat_results(p,b)

if nargin==0
    p = parameters();
    b = var_bounds();
end

X_nom = [b.X_noms; 1];
idx_F = strcmp(b.var_names, 'F_max');
F_nom = X_nom(idx_F);

force_lim_frac = linspace(0.1,1.1,30);
[power,force,runtime] = deal(zeros(size(force_lim_frac)));

for i=1:length(force_lim_frac)
    X = X_nom;
    X(idx_F) = force_lim_frac(i) * F_nom;
    tic
    [~, ~, ~, ~, val] = simulation(X,p);
    runtime(i) = toc;
    power(i) = val.power_avg;
    force(i) = val.force_heave_op;
end

fig = figure;
yyaxis left
plot(force_lim_frac * F_nom, power/1e3)
ylabel('Average Electrical Power (kW)')
hold on
yyaxis right
plot(force_lim_frac * F_nom, force/1e6)
plot([1 1]*F_nom, ylim, 'k--')
ylabel('Structural Fatigue Load (MN)')
xlabel('PTO Force Limit (MN)')
improvePlot

figure
plot(force_lim_frac * F_nom, runtime*1e3)
ylabel('Simulation Runtime (ms)')
xlabel('PTO Force Limit (MN)')
improvePlot

end

