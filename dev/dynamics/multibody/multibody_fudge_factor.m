p = parameters();
b = var_bounds();
X = [b.X_noms; 1];
wecsim_filename = 'wecsim_sparcd0_floatcd0_14e825b_01a25103-1490-4189-b376-76cbdf616634';
report = true;
override = false;
power_matrix_compare(X, p, wecsim_filename, report, override)

% add breakpoint on line 115 of power_matrix_compare and do the following in
% command window: (see 1/27/25 slides for screenshots)
sp = sim.power_mech_unsat;
ap = actual.power_mech_unsat;
ratio = ap ./ sp

figure; plot(T,ratio)
ylabel('RM3 report power / wecsim singlebody power')
xlabel('T')

y=ratio(H==5.75,:);
figure; plot(T,y)

% then do apps > curve fitter with rational function: num degree 0, den degree 2