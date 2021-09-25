function [F_ptrain,P_elec] = ptrain(s, u, i_PT)

% in reality these should be functions of i_PT
F_max = Inf;
eff = 1;

F_ptrain = min(u,F_max);
P_elec = F_ptrain * s(2) * eff;

end

