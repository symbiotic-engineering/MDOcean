function [feasible, A_ineq, failed, feasible_lin] = is_feasible(g_nonlin, x, p, b, idx_ignore)

if nargin<5
    idx_ignore = false(size(g_nonlin));
end

tol = -.01;
feasible_nonlin = all(g_nonlin(~idx_ignore) >= tol);

% The linear constraint currently assume a certain design variable order.
% If you change the order, this assert reminds you to update A_ineq & b_ineq.
assert(all( strcmp(b.var_names(1:4),{'D_s','D_f','T_f_2','h_s'}) ))

[A_ineq, b_ineq] = lin_ineq_constraints(p);
if isrow(x)
    x = x.';
end
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
        if g_linear(i) < 0
            failed = [failed ', ' lin_const_names{i}];
        end
    end

end

end