
function [F_max, D_env, P_elec] = dynamicSimulation(x,p,m_float)

% run simulation
sol = ode45(@(t,s)dynamics(t,s,x,p,m_float), [0 p.tfinal], p.s0);
time = 0:p.dt:p.tfinal;
state = deval(sol,time);
[~, P_elec, D_env, F] = dynamics(time, state, x, p, m_float);

% plot results
figure
plot(time,state*1e4,time,P_elec/100,time,D_env,time,F)
legend({'position*1e4','velocity*1e4','P elec/100','D env','net force'})
xlabel('time')

% covert time series to scalar outputs
D_env = mean(D_env); 
F_max = max(F);
P_elec = -mean(P_elec);

end

function [sdot, P_elec, D_env, F] = dynamics(t, s, x, p, m_float)
% s    = [position velocity]
% sdot = [velocity acceleration]

u = controls(s, x.D_int);
[F_ptrain, P_elec] = ptrain(s, u, p.i_PT);
[F_hydro, m, D_env] = hydro(t, s, m_float, x.D_sft, p.d_WEC, x.N_WEC, ...
                            p.d_farm, p.d_shore, p.rho_w, p.g, p.Hs, p.T);

F = F_hydro + F_ptrain;

sdot = zeros(2,size(s,2));
sdot(1,:) = s(2,:);
sdot(2,:) = F./m;

end

