function [X,vals] = compare_run(p,b)

if nargin==0
    p = parameters();
    b = var_bounds();
end

x0_input = b.X_start_struct;

[Xs_opt,~,~,~,~,~,~,vals_opt] = gradient_optim(x0_input,p,b);
X_minLCOE = Xs_opt(:,1);
X_minCapex = Xs_opt(:,2);
val_minLCOE = vals_opt(1);
val_minCapex = vals_opt(2);

[X_maxPower,val_maxPower] = max_avg_power(p,b,X_minLCOE);

p_bal = p;
p_bal.avg_power_min = b.power_balanced;
[X_balanced,~,~,~,~,~,~,val_balanced] = gradient_optim(x0_input,p_bal,b,2);

X_nom	= [b.X_noms' 1];
[~,~,~,val_nom] = simulation(X_nom,p);

%%
X    = [X_nom;   X_minLCOE';  X_minCapex';  X_maxPower';  X_balanced'];
vals = [val_nom, val_minLCOE, val_minCapex, val_maxPower, val_balanced];

end