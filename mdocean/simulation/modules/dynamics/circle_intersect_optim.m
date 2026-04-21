function [p_star,candidates,constraints_active] = circle_intersect_optim(c,r)
% finds the point p_star=[x,y] closest to the origin that is 
% within the intersection of circular discs with centers c and radii r.
% This is equivalent to the following optimization problem: 
% p_star = argmin_p |p| s.t. |p-c_i|-r_i <= 0 for all i.
% If the origin is feasible, p_star is the origin.
% If not, p_star will either be on the intersections of the circles,
% or on the point on each circle closest to the origin - these points
% are returned as candidates. If the problem is infeasible, 
% p_star is returned empty, constraints_active is all true, and a warning is issued.

N = length(r);
max_num_intersects = N * (N-1);
max_num_candidates = max_num_intersects + N;
candidates = NaN(max_num_candidates,2);
circle_indices = zeros(max_num_candidates,2);

% 1) check if origin feasible
origin = [0,0];
origin_feasible = check_feasible(origin,c,r);
if origin_feasible
    p_star = origin;
    constraints_active = false(size(r));
else

    % 2) pairwise intersections
    count = 0;
    for i = 1:N
        for j = i+1:N
            [xs, ys] = circcirc(c(i,1),c(i,2),r(i),c(j,1),c(j,2),r(j));
            circle_indices(count+(1:2),:) = [i,j;i,j];
            candidates(count+(1:2),:) = [xs.', ys.']; 
            count = count + 2;
        end
    end
    
    % 3) single-circle closest boundary
    for i = 1:N
        c_norm = norm(c(i,:));
        if c_norm == 0
            p = [r(i), 0];
        else
            p = c(i,:) - r(i) * (c(i,:) / c_norm);
        end
        candidates(count+i,:) = p;
        circle_indices(count+i,:) = [i,i];
    end
    
    % 4) filter feasible
    feasible_all_circles = false(length(candidates),1);
    for i=1:length(candidates)
        p = candidates(i,:);
        feasible_all_circles(i) = check_feasible(p,c,r);
    end
    feasible_p = candidates(feasible_all_circles,:);
    
    % 5) choose closest
    [~,idx] = min(vecnorm(feasible_p,2,2));
    p_star = feasible_p(idx,:);
    
    if isempty(p_star)
        constraints_active = true(size(r));
        if any(isnan(candidates))
            idx_nan = any(isnan(candidates),2);
            circles_no_intersect = unique(circle_indices(idx_nan,:),'rows');
            for i = 1:size(circles_no_intersect,1)
                warning('infeasible: circle %i and %i do not intersect',circles_no_intersect(i,1),circles_no_intersect(i,2))
            end
        else
            warning('each constraint pair is individually feasible, but all constraints together are infeasible')
        end
    else
        constraints_active = abs(eval_constraint(p_star,c,r)) < 1e-4;
    end
end
end

function feasible_all_circles = check_feasible(p,c,r)
% checks if a single point p is within all circles
% returns a scalar boolean
    dist_outside = eval_constraint(p,c,r);
    feas_each_circle = dist_outside <= 1e-4;
    feasible_all_circles = all(feas_each_circle);
end

function constr_fcn = eval_constraint(p,c,r)
% checks how far a single point p is in violation of the constraint
% positive means not ok (outside circle), negative means ok (inside circle)
% returns a scalar value
    assert(all(size(p)==[1 2]))
    constr_fcn = vecnorm(p - c,2,2) - r;
end
