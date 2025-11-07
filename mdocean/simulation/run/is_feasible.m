function [feasible, A_ineq, failed, feasible_lin] = is_feasible(g_nonlin, x, p, b, idx_ignore)

if nargin<5
    idx_ignore = false(size(g_nonlin));
end

tol = -.01;
feasible_nonlin = all(g_nonlin(~idx_ignore) >= tol);

% The linear constraint currently assume a certain design variable order.
% If you change the order, this assert reminds you to update A_ineq & b_ineq.
assert(all( strcmp(b.var_names,{'D_s','D_f','T_f_2','h_s','h_fs_clear',...
    'F_max','P_max','t_fb','t_sr','t_d','h_stiff_f','h1_stiff_d','M'}) ))

[A_ineq, b_ineq] = lin_ineq_constraints(p);
if isrow(x)
    x = x.';
end
g_linear = b_ineq-A_ineq*x;

feasible_lin = all(g_linear >= 0);

feasible = feasible_nonlin && feasible_lin;

if nargout > 2
    nl_const_names = b.constraint_names;
    lin_const_names = b.lin_constraint_names;

    failed_nl  = strjoin(nl_const_names(g_nonlin  < tol),', ');
    failed_lin = strjoin(lin_const_names(g_linear < tol),', ');
    if isempty(failed_lin) || isempty(failed_nl)
        delim = '';
    else
        delim = ', ';
    end
    failed = strjoin({failed_nl,failed_lin},delim);

end

end
