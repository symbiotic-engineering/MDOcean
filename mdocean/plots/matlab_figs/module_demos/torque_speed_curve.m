function [fig] = torque_speed_curve(X,p,b,color)

F_max = X(strcmp(b.var_names, 'F_max')) * 1e6;
P_max = X(strcmp(b.var_names, 'P_max')) * 1e5;

[~,~,~,val] = simulation(X, p);


V_u = val.X_u .* val.w;
V_max = max(V_u(:));

% points: v [m/s], F [N]
stall_point = [0, F_max];
base_point = [P_max/F_max, F_max];
max_point = [V_max, P_max/V_max];


end