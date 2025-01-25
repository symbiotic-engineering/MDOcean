function [feasible, failed] = is_feasible(g_nonlin, x, p, b)

tol = -.01;
feasible_nonlin = all(g_nonlin >= tol);

% The linear constraint currently assume a certain design variable order.
% If you change the order, this assert reminds you to update A_ineq & b_ineq.
assert(all( strcmp(b.var_names(1:4),{'D_s','D_f','T_f_2','h_s'}) ))

A_ineq = [-p.D_d_over_D_s 0 0                        0;
            1            -1 0                        0;
          -p.T_s_over_D_s 0 1                        0
          p.T_s_over_D_s, 0, 1/p.T_f_2_over_h_f - 1, -1
            0             p.harmonics*pi/(2*700.5) 1 0
            p.harmonics*pi/(2*700.5) + p.T_s_over_D_s, 0 0 0];
b_ineq = [-p.D_d_min -0.01 -0.01 -0.01 p.h p.h]';
g_linear = b_ineq-A_ineq*x(1:4);

feasible_lin = all(g_linear >= 0);

feasible = feasible_nonlin && feasible_lin;

if nargout > 1
    nl_const_names = b.constraint_names;
    lin_const_names = b.lin_constraint_names;

    failed = '';
    for i = 1:length(g_nonlin)
        if g_nonlin(i) < tol
            failed = [failed ', ' nl_const_names{i}];
        end
    end
    for i = 1:length(g_linear)
        if g_linear < 0
            failed = [failed ', ' lin_const_names{i}];
        end
    end

end

end