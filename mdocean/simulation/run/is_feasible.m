function [feasible, A_ineq, failed] = is_feasible(g_nonlin, x, p, b)

tol = -.01;
feasible_nonlin = all(g_nonlin >= tol);

A_ineq = [-p.D_d_over_D_s 0 0                        0;
            1            -1 0                        0; % fixme use D_f_in
          -p.T_s_over_D_s 0 1                        0
          p.T_s_over_D_s, 0, 1/p.T_f_2_over_h_f - 1, -1
            0             p.harmonics*pi/(2*700.5) 1 0
            p.harmonics*pi/(2*700.5) + p.T_s_over_D_s, 0 0 0];
b_ineq = [-p.D_d_min -0.01 -0.01 -0.01 p.h p.h]';
g_linear = b_ineq-A_ineq*x(1:4);

feasible_lin = all(g_linear >= 0);

feasible = feasible_nonlin && feasible_lin;

if nargout > 2
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