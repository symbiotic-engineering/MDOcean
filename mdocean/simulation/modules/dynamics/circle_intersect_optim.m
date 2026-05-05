function [p_star,candidates,constraints_active] = circle_intersect_optim(c,r,warn_if_infeasible)
% Finds the point p_star=[x,y] closest to the origin that is 
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
too_far_check = false(max_num_candidates,1);

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
            [xs, ys, too_far] = circle_circle_intersect(c(i,1),c(i,2),r(i),c(j,1),c(j,2),r(j));
            circle_indices(count+(1:2),:) = [i,j;i,j];
            too_far_check(count+(1:2)) = too_far;
            candidates(count+(1:2),:) = [xs.', ys.']; 
            count = count + 2;
        end
    end
    
    % 3) single-circle closest boundary
    for i = 1:N
        dist = norm(c(i,:));
        if dist < 1e-10
            % center at origin: any boundary point is equally close
            p = [r(i), 0];
        else
            p = c(i,:) - r(i) * (c(i,:) / dist);
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
        if any(isnan(candidates(:)))
            idx_nan = any(isnan(candidates),2);
            idx_infeasible = idx_nan & too_far_check; % this assumes all circles are feasible inside, infeasible outside
            circles_no_intersect = unique(circle_indices(idx_infeasible,:),'rows');
            if warn_if_infeasible
                for i = 1:size(circles_no_intersect,1)
                    warning('MDOcean:circle_intersect_optim:infeasible',...
                        'infeasible: circle %i and %i do not intersect and are not contained/coincident',...
                        circles_no_intersect(i,1),circles_no_intersect(i,2))
                end
            end
        else
            if warn_if_infeasible
                warning('MDOcean:circle_intersect_optim:infeasible',...
                    ['Each constraint pair is individually feasible, '...
                    'but all constraints together are infeasible'])
            end
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

function [xs, ys, too_far] = circle_circle_intersect(x1,y1,r1,x2,y2,r2)
% Finds the intersection points of two circles.
% Returns NaN if circles don't intersect or are coincident.
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
    
    too_far = d > r1 + r2;
    contained_or_coincident = d < abs(r1-r2) || d == 0;
    if too_far || contained_or_coincident
        % no intersection (too far, contained, or coincident)
        xs = [NaN, NaN];
        ys = [NaN, NaN];
        return
    end
    
    a = (r1^2 - r2^2 + d^2) / (2*d);
    h = sqrt(max(r1^2 - a^2, 0));  % max with 0 for numerical safety
    
    % midpoint
    mx = x1 + a*(x2-x1)/d;
    my = y1 + a*(y2-y1)/d;
    
    % offset
    dx = h*(y2-y1)/d;
    dy = h*(x2-x1)/d;
    
    xs = [mx+dx, mx-dx];
    ys = [my-dy, my+dy];
end
