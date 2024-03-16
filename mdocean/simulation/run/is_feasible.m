function [feasible, failed] = is_feasible(g, b)

tol = -.01;
feasible = all(g >= tol);

if nargout > 1
    const_names = b.constraint_names;

    % delete once I find where b.constraint_names is defined and add more
    % names
    const_names_temp = const_names
    if length(const_names) < length(g)
        for i = length(const_names)+1:length(g)
            const_names_temp{i} = 'prevent slamming';
        end
    end
    const_names = const_names_temp;
    failed = ' ';
    for i = 1:length(g)
        if g(i) < tol
            failed = [failed const_names{i}];
        end
    end
end

end