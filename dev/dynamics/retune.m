% retune heave force and power multiplier to RM3 report
p = parameters('report');
b = var_bounds('report');
[~,val] = run_single(p,b);
actual = validation_inputs('report');

current_vals = [   val.force_heave,    val.power_avg];
desired_vals = [actual.force_heave, actual.power_avg];
current_mult = [p.F_heave_mult p.power_scale_multibody]
desired_mult = current_mult .* desired_vals ./ current_vals