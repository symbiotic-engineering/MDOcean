function [feasible, failed] = is_feasible(g, b)

tol = -.01;
feasible = all(g >= tol);

if nargout > 1
    const_names = b.constraint_names;

    failed = ' ';
    for i = 1:length(g)
        if g(i) < tol
            failed = [failed const_names{i}];
        end
    end
end

end