function [F_ptrain,P_elec] = ptrain(s, u, i_PT)

F_max = Inf; % in reality this is a function of i_PT
F_ptrain = min(u,F_max);

eff = 1; % in reality this is a function of i_PT, s, and F_ptrain
P_elec = F_ptrain .* s(2,:) * eff;

end

