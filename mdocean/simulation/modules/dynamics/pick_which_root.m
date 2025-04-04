function mult = pick_which_root(roots, idx_no_sat, idx_zero, a_quad, b_quad, c_quad)

    which_soln = roots == real(roots) & roots > 0 & roots <= 1; % real solns on (0, 1]

    both_ok = sum(which_soln,3) == 2;
    
    which_soln(idx_no_sat | idx_zero) = 1;  % temporarily mark the limit solutions
                                            % as having one solution, to ensure the 
                                            % logic below works correctly

    if any(both_ok,'all')       % two solutions
        roots_equal = roots(:,:,1) == roots(:,:,2);
        if ~all(roots_equal(both_ok),'all')
            warning('Two valid non-equal solutions, choosing 1st arbitrarily')
        end

        mult = handle_two_solns(both_ok,which_soln,roots,idx_no_sat,idx_zero,a_quad,b_quad,c_quad);
    
    else                        % one or zero solutions
        num_solns = sum(which_soln,3);
        if ~(all( num_solns == 1,'all') ) % confirm that 1 soln per sea state meets criteria
            % if there are some with no solns, proceed with which_soln set to zero
            % for the problem sea states, which will set mult = 0
            warning('Some sea states have no valid quadratic solution, so their energy is zeroed.')        
        end

        mult = get_relevant_soln(which_soln,roots,idx_no_sat, idx_zero);   
    end
end

function mult = get_relevant_soln(which_soln, roots, idx_no_sat, idx_zero)
% pick the specified roots using multidimensional logical indexing
    mult = zeros(size(idx_no_sat));

    % figure out 3d and 2d indices
    idx_3d_first_sol = which_soln;
    idx_3d_first_sol(:,:,2) = false;
    idx_3d_second_sol = which_soln;
    idx_3d_second_sol(:,:,1) = false;
    idx_2d_first_sol = which_soln(:,:,1);
    idx_2d_second_sol = which_soln(:,:,2);
    
    mult(idx_2d_first_sol) = roots(idx_3d_first_sol);
    mult(idx_2d_second_sol) = roots(idx_3d_second_sol);
    mult(idx_no_sat) = 1;
    mult(idx_zero) = 0;
end

function mult = handle_two_solns(both_ok, which_soln, roots, idx_no_sat,idx_zero,a,b,c)
     
%     % analytical criteria for root between 0 and 1
%     use_1 = (a >= 0 & c <= 0 & (a+b+c) <= 0) | ...
%             (a <= 0 & c <= 0 & (a+b+c) >= 0);
%     use_2 = (a <= 0 & c >= 0 & (a+b+c) <= 0)| ...
%             (a >= 0 & c == 0 & (a+b+c) >= 0);


    [row,col] = find(both_ok);

    % if both are ok, choose the first one arbitrarily for now
    which_soln(row,col,2) = false; 
    mult_1 = get_relevant_soln(which_soln,roots,idx_no_sat, idx_zero);   

    mult = mult_1;
% 
%     % now choose the second one arbitrarily
%     which_soln_2 = which_soln;
%     which_soln_2(row,col,2) = true;
%     which_soln_2(row,col,1) = false;
%     mult_2 = get_relevant_soln(which_soln_2,roots,idx_no_sat);   
% 
%     % compare outliers to figure out if first or second is better
%     [use_1_vals,use_2_vals] = compare_outliers(mult_1,mult_2,both_ok);
%     [use_1_idxs,use_2_idxs] = compare_outliers(double(which_soln(:,:,1)),...
%                                             double(which_soln_2(:,:,1)),both_ok);    
%     
%     use_1 =  (use_1_vals || use_1_idxs) && ~(use_2_vals || use_2_idxs);
%     use_2 = ~(use_1_vals || use_1_idxs) &&  (use_2_vals || use_2_idxs);
% 
%     if use_1
%         mult = mult_1;
%     elseif use_2
%         mult = mult_2;
%     else
%         figure
%         subplot 121
%         contourf(mult_1)
%         subplot 122
%         contourf(mult_2)
% 
%         error(['Failed to figure out which solution to the quadratic ' ...
%             'equation is relevant, try manual inspection.'])
%     end
end

function [use_1,use_2] = compare_outliers(array_1,array_2,relevant_idx)
    window = 6;
    outliers_1 = isoutlier(array_1,'movmedian',window);
    outliers_2 = isoutlier(array_2,'movmedian',window);
    
    outliers_1_relevant = outliers_1(relevant_idx);
    outliers_2_relevant = outliers_2(relevant_idx);

    one_all_outliers = all(outliers_1_relevant);
    two_all_outliers = all(outliers_2_relevant);
    one_all_ok = all(~outliers_1_relevant);
    two_all_ok = all(~outliers_2_relevant);
    
    use_1 = (two_all_outliers && ~one_all_outliers) || (one_all_ok && ~two_all_ok);
    use_2 = (one_all_outliers && ~two_all_outliers) || (two_all_ok && ~one_all_ok);
%     use_1 = sum(outliers_1_relevant,'all') < sum(outliers_2_relevant,'all');
%     use_2 = sum(outliers_1_relevant,'all') > sum(outliers_2_relevant,'all');
end

