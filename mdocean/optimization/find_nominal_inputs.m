% find the values of B_p, w_n, and F_max for the nominal RM3 
% that make P_elec and F_hydro match the report values

% setup
clear;clc;close all
p = parameters();
p.N_WEC = 1;
b = var_bounds(p);
v = validation_inputs();

p.power_max = v.power_max;

y_desired = [v.power_avg, v.force_heave, 1];

% x = [F_max_nom, B_p_nom, w_n_nom]
x_min = [b.F_max_min b.B_p_min b.w_n_min];
x_max = [b.F_max_max b.B_p_max b.w_n_max];
x0 = [5 10 2*pi/8];
fcn = @(x) errFunc(x,y_desired,p,b);

% optimize
x = fmincon(fcn, x0,[],[],[],[],x_min,x_max);

% check feasibility
X = [b.X_noms; 1];
X(5:7) = x;
[LCOE, P_var, ~, g] = simulation(X,p);
[feasible,failed] = is_feasible(g, b)

% display x output
array2table(x,'VariableNames',{'F_max (1e6 N)','B_p (1e6 Ns/m)','w_n (rad/s)'})

% display y output
[~,y] = errFunc(x,y_desired,p,b);
results = round([y; y_desired] ./ [1000 1000 1],3,'significant');
array2table(results,'RowNames',{'Sim Output','RM3 Actual'},...
    'VariableNames',{'Average Power (kW)','Max Structural Force (kN)','Max Powertrain Force Ratio (-)'})

function [err,y] = errFunc(x,y_desired,p,b)
    X = [b.X_noms; 1];
    X(5:7) = x;
    [~, ~, ~, g, val] = simulation(X, p);
    y = [val.power_avg, val.force_heave, g(14)+1];
    err = abs(y - y_desired) ./ y_desired;
    err = norm(err);
end
    
