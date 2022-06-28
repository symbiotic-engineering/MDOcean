% find the values of B_p, w_n, and F_max for the nominal RM3 
% that make P_elec, max(P_matrix), and F_hydro match the report values

clear;clc;close all
p = parameters();
b = var_bounds(p);

P_elec_desired = 286e3 * 0.3;
P_max_desired = 286e3;
F_hydro_desired = 8500e6;
y_desired = [P_elec_desired, P_max_desired, F_hydro_desired];

% x = [F_max_nom, B_p_nom, w_n_nom]
x_min = [b.F_max_min b.B_p_min b.w_n_min];
x_max = [b.F_max_max b.B_p_max b.w_n_max];
x0 = [5 10 2*pi/8];
fcn = @(x) errFunc(x,y_desired,p,b);
x = fmincon(fcn, x0,[],[],[],[],x_min,x_max);

errFunc(x,y_desired,p,b)

function err = errFunc(x,y_desired,p,b)
    X = [b.X_noms; 1];
    X(5:7) = x;
    [~, ~, ~, ~, ~, ~, ~, ~, ...
        P_elec, ~, P_matrix, ~, val] = simulation(X, p);
    y = [P_elec, max(P_matrix,[],'all'), val(7)];
    err = abs(y - y_desired) ./ y_desired;
    err = norm(err);
end
    
