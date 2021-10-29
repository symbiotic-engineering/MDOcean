
function [F_heave, F_surge, F_ptrain, D_env, P_elec] = dynamicSimulation(x,p,m_float,t_f)

% run simulation
sol = ode45(@(t,s)dynamics(t,s,x,p,m_float,t_f), [0 p.tfinal], p.s0);
time = 0:p.dt:p.tfinal;
state = deval(sol,time);
[~, P_elec, D_env, F_heave, F_surge, F_ptrain] = dynamics(time, state, x, p, m_float, t_f);

% plot results
% figure
% plot(time,state*1e4,time,P_elec/100,time,D_env,time,F)
% legend({'position*1e4','velocity*1e4','P elec/100','D env','net force'})
% xlabel('time')

% covert time series to scalar outputs
D_env = mean(D_env); 
F_heave = max(F_heave);
F_surge = max(F_surge);
F_ptrain = max(F_ptrain);
P_elec = -mean(P_elec);

end

function [sdot, P_elec, D_env, F_heave, F_surge, F_ptrain] = dynamics(t, s, x, p, m_float, t_f)
% s    = [position velocity]
% sdot = [velocity acceleration]

u = controls(s, x.D_int);
[F_ptrain, P_elec] = ptrain(s, u, p.i_PT);
[F_heave, F_surge, m, D_env] = hydro(t, s, m_float, x.D_sft,x.D_i,x.D_or,p.t_r, p.rho_w, p.g, p.Hs, p.T, t_f);
F_net = F_heave + F_ptrain;

sdot = zeros(2,size(s,2));
sdot(1,:) = s(2,:);
sdot(2,:) = F_net./m;

end

