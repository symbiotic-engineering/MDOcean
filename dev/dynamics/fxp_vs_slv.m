fxp_iters = [10:10:90];
slv_iters = [10:10:90];

p = parameters();
%b = var_bounds();
X = [b.X_noms; 1];

[FXP,SLV] = meshgrid(fxp_iters,slv_iters);
time = zeros(size(FXP));
for i=1:numel(FXP)
    p.fxp_max_iters = FXP(i);
    p.slv_max_iters = SLV(i);
    t = tic;
    [~,~,sol] = simulation(X, p);
    time(i) = toc(t);
end

f = figure;
contourf(FXP, SLV, time, 20);
colorbar
xlabel('Fixed-Point Iterations');
ylabel('Solver Iterations');
title('Simulation Time (s)');
save_pdf(f, 'fxp_vs_slv.pdf');