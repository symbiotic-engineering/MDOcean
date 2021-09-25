
function [F_max, D_env, P_elec] = dynamicSimulation(x,p)

% run simulation
sol = ode45(@(t,s)dynamics(t,s,x,p), [0 p.tfinal], p.s0);
time = 0:p.dt:p.tfinal;
state = deval(sol,time);
%[~, P_elec, D_env, F] = dynamics(0, state, x, p);

% plot results
figure
plot(time,state)%time,P_elec,time,D_env,time,F)
legend({'position','velocity'})%'P elec','D env','F')

% covert time series to scalar outputs
D_env = 0;%mean(D_env); 
F_max = 0;%max(F);
P_elec = 0;%mean(P_elec);

end

function [sdot, P_elec, D_env, F] = dynamics(t, s, x, p)
% s    = [position velocity]
% sdot = [velocity acceleration]

u = controls(s, x.D_int);
[F_ptrain, P_elec] = ptrain(s, u, x.i_PT);
[F_hydro, m, D_env] = hydro(t, s, x.L, x.W, x.M, x.d_WEC, x.N_WEC, ...
                            x.d_farm, p.d_shore, p.rho_w, p.g, p.Hs, p.T);

F = F_hydro + F_ptrain;

sdot = [0; 0];
sdot(1) = s(2);
sdot(2) = F/m;

end

