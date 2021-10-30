function [F_heave, F_surge, F_ptrain, D_env, P_elec] = dynamicSimulation(x,p,m_float,t_f)

time = 0:p.dt:p.tfinal;

num_Hs = length(p.Hs);
num_T = length(p.T);
P_elec = zeros(num_Hs,num_T);

for H_idx = 1:num_Hs
    for T_idx = 1:num_T
        if p.JPD(H_idx,T_idx) ~= 0
            Hs = p.Hs(H_idx);
            T = p.T(T_idx);

            % run simulation
            sol = ode45(@(t,s)dynamics(t,s,x,p,Hs,T,m_float,t_f), [0 p.tfinal], p.s0);
            state = deval(sol,time);
            [~, P_elec_vs_time, D_env, F_heave, F_surge, F_ptrain] = dynamics(time, state, x, p, Hs, T, m_float, t_f);
            P_elec(H_idx,T_idx) = -mean(P_elec_vs_time);

            sol = ode45(@(t,s)dynamics(t,s,x,p,Hs,T,m_float,t_f), [0 p.tfinal], p.s0);
            state = deval(sol,time);
            [~, ~, ~, F_heave, F_surge, F_ptrain] = dynamics(time, state, x, p, p.Hs_struct, p.T_struct, m_float, t_f);
            
            % plot results
            % figure
            % plot(time,state*1e4,time,P_elec/100,time,D_env,time,F_heave,F_ptrain)
            % legend({'position*1e4','velocity*1e4','P elec/100','D env','net force'})
            % xlabel('time')
        end
    end 
end

% weight power across all sea states
P_weighted = P_elec .* p.JPD;
P_elec = mean(P_weighted,'all');

% covert time series to scalar outputs
D_env = mean(D_env); 
F_heave = max(F_heave);
F_surge = max(F_surge);
F_ptrain = max(F_ptrain);

end

function [sdot, P_elec, D_env, F_heave, F_surge, F_ptrain] = dynamics(t, s, x, p, Hs, T, m_float, t_f)
% s    = [position velocity]
% sdot = [velocity acceleration]

u = controls(s, x.D_int);
[F_ptrain, P_elec] = ptrain(s, u);
[F_heave, F_surge, m, D_env] = hydro(t, s, m_float, x.D_sft, p.rho_w, p.g, Hs, T, t_f);
F_net = F_heave + F_ptrain;

sdot = zeros(2,size(s,2));
sdot(1,:) = s(2,:);
sdot(2,:) = F_net./m;

end

