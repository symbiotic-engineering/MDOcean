% retune heave force and power multiplier to RM3 report
p = parameters('report');
p.control_type = 'damping';
b = var_bounds('report');
b.F_max_nom = find_nominal_inputs(b,p);
b.X_noms(strcmp(b.var_names,'F_max')) = b.F_max_nom;

[~,~,val] = run_single(p,b);
actual = validation_inputs('report');

current_vals = [   val.force_heave,    val.power_avg];
desired_vals = [actual.force_heave, actual.power_avg];
current_mult = [p.F_heave_mult p.power_scale_multibody]
desired_mult = current_mult .* desired_vals ./ current_vals