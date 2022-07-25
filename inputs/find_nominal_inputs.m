% find the values of B_p, w_n, and F_max for the nominal RM3 
% that make P_elec, max(P_matrix), and F_hydro match the report values

% setup
clear;clc;close all
p = parameters();
b = var_bounds(p);
v = validation_inputs();

p.power_max = v.power_max;

y_desired = [v.power_avg, v.power_max, v.force_heave];

% x = [F_max_nom, B_p_nom, w_n_nom]
x_min = [b.F_max_min b.B_p_min b.w_n_min];
x_max = [b.F_max_max b.B_p_max b.w_n_max];
x0 = [5 10 2*pi/8];
fcn = @(x) errFunc(x,y_desired,p,b);

% optimize
x = fmincon(fcn, x0,[],[],[],[],x_min,x_max)

% check feasibility
X = [b.X_noms; 1];
X(5:7) = x;
[LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM, P_elec, D_d, ~, g] = simulation(X,p);
FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
[feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p)

% display outputs
[~,y] = errFunc(x,y_desired,p,b);
results = round([y; y_desired] / 1000);
tab = array2table(results,'RowNames',{'Sim Output','RM3 Actual'},...
    'VariableNames',{'Average Power (kW)','Max Power (kW)','Max Force (kN)'})

function [err,y] = errFunc(x,y_desired,p,b)
    X = [b.X_noms; 1];
    X(5:7) = x;
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, val] = simulation(X, p);
    y = [val.power_avg, val.power_max, val.force_heave];
    err = abs(y - y_desired) ./ y_desired;
    %err(1:2) = 0; % override max power
    err = norm(err);
end
    
