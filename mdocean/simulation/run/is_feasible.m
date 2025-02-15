function [feasible, A_ineq, failed] = is_feasible(g_nonlin, x, p, b)

tol = -.01;
feasible_nonlin = all(g_nonlin >= tol);

[A_ineq, b_ineq] = lin_ineq_constraints(p);
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