


function [sigma_vm_bot,sigma_vm_top] = float_plate_stress(D_f, D_f_in, F_heave, num_sections, t_bot, t_top, h_stiff, w_stiff, D_f_tu, nu)

    A = pi/4 * (D_f^2 - D_f_in^2);
    P_hydrodynamic = F_heave / A;
    W = F_heave / num_sections;
    
    b_out = pi * D_f    / num_sections;
    b_in  = pi * D_f_in / num_sections;
    b_avg = sqrt(b_out * b_in);
    
    h = (D_f-D_f_in)/2;

    % trapezoid interpolation
    % see page 29 of notebook 12/15/24
    m = (b_out-b_in)/(2*h);
    if h >= b_avg
        shorter_length = b_in*(m+sqrt(1+m^2));
        longer_length = h;
    else
        shorter_length = h;
        longer_length  = (b_out+b_in)/2 + h*(1-sqrt(1+m^2));
    end
    
    length_ratio = longer_length / shorter_length;
    r0 = 1.02;%D_f_tu / 2;
    width_plate = b_out; % fixme, should be combo of b_out and b_in? but b_in gives imag number

    sigma_vm_bot = bottom_plate_stress(length_ratio, shorter_length, P_hydrodynamic, t_bot, h_stiff, w_stiff, width_plate);
    sigma_vm_top = top_plate_stress(length_ratio, shorter_length, W, t_top, h_stiff, w_stiff, r0, width_plate, nu);

end

function sigma_vm_top = top_plate_stress(length_ratio, shorter_length, W, t_top, h_stiff, w_stiff, r0, width_plate, nu)

    length_ratio_data = [1 : 0.2 : 2 1000];
    beta_1_data = [-0.238 -0.078 0.011 0.053 0.068 0.067 0.067];
    beta_2_data = [0.7542 0.8940 0.9624 0.9906 1.0000 1.004 1.008];
    beta_1_fcn  = @(length_ratio) interp1(length_ratio_data, beta_1_data, length_ratio);
    beta_2_fcn  = @(length_ratio) interp1(length_ratio_data, beta_2_data,  length_ratio);
    beta_1 = beta_1_fcn(length_ratio);
    beta_2 = beta_2_fcn(length_ratio);

    % Roark's table 11.4 case 8b page 514
    % fixme - this is potentially not applicable since it assumes force in
    % middle, but my force is on the edge and therefore creates ~no bending moment
    % perhaps check ref 26 and 21 to see if they have a more general (off center) equation 
    sigma_edge_no_stiff = 3*W/(2*pi*t_top^2) * ( (1+nu)*log(2*shorter_length/(pi*r0)) + beta_1);
    sigma_cent_no_stiff = -beta_2 * W / t_top^2;

    M_edge = sigma_edge_no_stiff * t_top^2 / 6;
    M_cent = sigma_cent_no_stiff * t_top^2 / 6;
    
    M_max = max(abs([M_edge,M_cent]));
    sigma_max = stiffened_plate_stress(t_top, h_stiff, width_plate, w_stiff, M_max);
    sigma_vm_top = sigma_max;

end

function sigma_vm = bottom_plate_stress(length_ratio, shorter_length, P_hydrodynamic, t_bot, h_stiff, w_stiff, width_plate)

    % data from table: Timoshenko table 35 p219
    length_ratio_data = [1 : 0.1 : 2 realmax];
    if length_ratio > length_ratio_data(end)
        error(['float structures: plate aspect ratio ' num2str(length_ratio) ' is too high.'])
    end
    alpha_shorter_data = -0.0001*[513 581 639 687 726 757 780 799 812 822 829 833];
    alpha_longer_data  = -0.0001*[513 538 554 563 568 570 571 571 571 571 571 571];
    alpha_shorter_fcn = @(length_ratio) interp1(length_ratio_data, alpha_shorter_data, length_ratio);
    alpha_longer_fcn  = @(length_ratio) interp1(length_ratio_data, alpha_longer_data,  length_ratio);
    
    M_shorter = alpha_shorter_fcn(length_ratio) * P_hydrodynamic * shorter_length^2;
    M_longer  = alpha_longer_fcn(length_ratio)  * P_hydrodynamic * shorter_length^2;
    
    sigma_shorter = stiffened_plate_stress(t_bot, h_stiff, width_plate, w_stiff, M_shorter);
    sigma_longer  = stiffened_plate_stress(t_bot, h_stiff, width_plate, w_stiff, M_longer);
    
    sigma_zz = P_hydrodynamic;
    
    % von mises
    s_11 = sigma_shorter;
    s_22 = 0;
    s_33 = sigma_zz;
    sigma_vm = sqrt(1/2 * ( (s_11 - s_22).^2 + (s_22 - s_33).^2 + (s_33 - s_11).^2 ));

end

function sigma = stiffened_plate_stress(t_plate, h_stiff, width_plate, width_stiff, moment_per_length)
    [h_eq,y_max] = get_stiffener_equivalent_properties(t_plate, h_stiff, width_plate, width_stiff);
    sigma = get_plate_stress(moment_per_length, y_max, h_eq);
end
