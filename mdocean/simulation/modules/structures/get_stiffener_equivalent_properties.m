function [h_eq,y_max,S_eq] = get_stiffener_equivalent_properties(t_plate, h_stiff, width_plate, width_stiff)
% https://ocw.mit.edu/courses/2-080j-structural-mechanics-fall-2013/f8fd2ad49d100766335b4e129a8a4791_MIT2_080JF13_Lecture7.pdf
% assumes stiffeners with rectangle shape (not I-beam)

    % renaming inputs to match MIT 2.080 equations
    h = t_plate;
    H = h_stiff;
    a = width_plate/2;
    b = width_stiff;
    
    % eq 7.67 from MIT 2.080, but  fixing typo: b/a should be b/2a
    % finding neutral axis
    num = 1 - b./(2*a) * (H/h)^2;
    den = 1 + b./(2*a) * (H/h);
    eta_over_h = 1/2 * num ./ den; 

    % eq 7.69 from MIT 2.080, but fixing typo: ^2 should be ^3 in term2, 
    % and eta should be eta over h in term1
    % finding equivalent height that creates equal area moment of inertia
    term1 = 1 - 3 * eta_over_h + 3 * eta_over_h.^2;
    term2 = (H/h)^3 + 3*(H/h)^2 * eta_over_h + 3 * H/h * eta_over_h.^2;
    h_eq_over_h_3 = 4 * ( term1 + b./(2*a) .* term2 ); 
    
    if any(h_eq_over_h_3 < 1) % if equivalent stiffened height is less than unstiffened height
        warning('stiffener equivalence broke')
        h_eq_over_h_3(h_eq_over_h_3 < 1) = 1;
    end
    h_eq_over_h = h_eq_over_h_3.^(1/3);
    h_eq = h_eq_over_h * h;
    
    if nargout > 1
        eta = eta_over_h * h; % location of neutral axis, from stiffened face of plate
        y_plate = h - eta; % neutral axis from unstiffened face of plate
        y_stiff = H + eta; % neutral axis from free face of stiffener
        y_max = max(y_plate, y_stiff); % max distance from neutral axis

        if nargout > 2   
            I = 1/12 * width_plate * h_eq.^3; % moment of inertia, of both the 
                                             % equivalent plate and the original stiffened plate
        
            S_eq = I ./ y_max; % equivalent section modulus
        end
    end
end