function sigma = get_stiffener_equivalent_stress(t_plate, h_stiff, width_plate, width_stiff, moment_per_length)
% https://ocw.mit.edu/courses/2-080j-structural-mechanics-fall-2013/f8fd2ad49d100766335b4e129a8a4791_MIT2_080JF13_Lecture7.pdf

    h = t_plate;
    H = h_stiff;
    a = width_plate/2;
    b = width_stiff;
    M = moment_per_length;
    
    eta_over_h = 1/2 * (1 - b/a * (H/h)^2) / (1 + b/a * (H/h)); % eq 7.67 MIT 2.080
    h_eq_over_h_3 = 4 * ( (1-3*eta_over_h+3*eta_over_h^2) + b/(2*a) * ( (H/h)^2 + 3*(H/h)^2*eta_over_h + 3*H/h*eta_over_h^2 )  ); % eq 7.69
    if h_eq_over_h_3 < 0
        warning('stiffener equivalence broke')
        h_eq_over_h_3 = 1;
    end
    h_eq_over_h = h_eq_over_h_3^(1/3);
    h_eq = h_eq_over_h * h;
    
    y_max = h-eta_over_h*h;

    sigma = 12 * M * y_max / h_eq^3;

end