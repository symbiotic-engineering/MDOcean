function [h_eq,y_max,S_eq] = get_stiffener_equivalent_properties(t_plate, h_stiff, width_plate, width_stiff)

    if     length(h_stiff)==1 && length(width_stiff)==1
        shape = 'tee';
    elseif length(h_stiff)==2 && length(width_stiff)==2
        shape = 'I';
    elseif length(h_stiff)==4 && length(width_stiff)==4
        shape = 'doubleI';
    else
        error('invalid input length')
    end

    % fixme if the stiffeners are spaced far apart, the effective width will be reduced due to shear lag
    if strcmp(shape,'tee')
        [centroid, h_eq_over_h_3, height_max] = T_beam_properties(t_plate, h_stiff, width_plate, width_stiff);

    elseif strcmp(shape,'I')
        [centroid, h_eq_over_h_3, height_max] = I_beam_properties(t_plate, h_stiff, width_plate, width_stiff);

    elseif strcmp(shape,'doubleI')
        [centroid, h_eq_over_h_3, height_max] = double_I_beam_properties(t_plate, h_stiff, width_plate, width_stiff);

    end

    if any(h_eq_over_h_3 < 1) % if equivalent stiffened height is less than unstiffened height
        warning('stiffener equivalence broke')
        h_eq_over_h_3(h_eq_over_h_3 < 1) = 1;
    end
    h = t_plate;
    h_eq_over_h = h_eq_over_h_3.^(1/3);
    h_eq = h_eq_over_h * h;

    if nargout > 1

        y_max = max(centroid, height_max - centroid); % max distance from neutral axis

        if nargout > 2
            I = 1/12 * width_plate * h_eq.^3; % moment of inertia, of both the
                                             % equivalent plate and the original stiffened plate

            S_eq = I ./ y_max; % equivalent section modulus
        end
    end
end

function [centroid, h_eq_over_h_3, height_max] = double_I_beam_properties(t_plate, h_stiff, width_plate, width_stiff)

    [M0,A0] = first_moment_of_area_rectangle(t_plate, width_plate, t_plate/2);
    [M1,A1] = first_moment_of_area_Ibeam(t_plate, h_stiff(1:2), width_stiff(1:2));
    [M2,A2] = first_moment_of_area_Ibeam(t_plate, h_stiff(3:4), width_stiff(3:4));

    moment = M0 + M1 + M2;
    area   = A0 + A1 + A2;
    centroid = moment./area; % centroid from bottom of plate

    I0 = second_moment_of_area_rectangle(t_plate, width_plate, 0);
    I1 = second_moment_of_area_Ibeam(t_plate, h_stiff(1:2), width_stiff(1:2));
    I2 = second_moment_of_area_Ibeam(t_plate, h_stiff(3:4), width_stiff(3:4));
    I = I0 + I1 + I2; % 2nd moment of area about bottom of plate

    h_eq_over_h_3 = equiv_I_beam(area, centroid, I, width_plate, t_plate);

    height_max = t_plate + max( sum(h_stiff(1:2)), sum(h_stiff(3:4)) );
end

function h_eq_over_h_3 = equiv_I_beam(area, centroid, I, width_plate, t_plate)
    % parallel axis theorem
    I_centroid = I - area .* centroid.^2; % 2nd moment of area about centroid

    h_eq_3 = 12 * I_centroid ./ width_plate;
    h_eq_over_h_3 = h_eq_3 / t_plate^3;

end

function I = second_moment_of_area_rectangle(height,width,height_start)
    I = width/3 * ( (height_start+height)^3 - (height_start)^3 );
end

function [centroid, h_eq_over_h_3, height_max] = I_beam_properties(t_plate, h_stiff, width_plate, width_stiff)
    [M0,A0] = first_moment_of_area_rectangle(t_plate, width_plate, t_plate/2);
    [M1,A1] = first_moment_of_area_Ibeam(t_plate, h_stiff(1:2), width_stiff(1:2));

    moment = M0 + M1;
    area   = A0 + A1;
    centroid = moment/area; % centroid from bottom of plate

    I0 = second_moment_of_area_rectangle(t_plate, width_plate, 0);
    I1 = second_moment_of_area_Ibeam(t_plate, h_stiff(1:2), width_stiff(1:2));
    I = I0 + I1; % 2nd moment of area about bottom of plate

    h_eq_over_h_3 = equiv_I_beam(area, centroid, I, width_plate, t_plate);

    height_max = t_plate + sum(h_stiff);
end

function I = second_moment_of_area_Ibeam(t_plate, h_stiff, width_stiff)
    I1 = second_moment_of_area_rectangle(h_stiff(1), width_stiff(1), t_plate);
    I2 = second_moment_of_area_rectangle(h_stiff(2), width_stiff(2), t_plate + h_stiff(1));
    I = I1 + I2;
end

function [M,A] = first_moment_of_area_Ibeam(t_plate, h_stiff, width_stiff)
    [M1,A1] = first_moment_of_area_rectangle(h_stiff(1), width_stiff(1), t_plate + 1/2*h_stiff(1));
    [M2,A2] = first_moment_of_area_rectangle(h_stiff(2), width_stiff(2), t_plate + h_stiff(1) + 1/2*h_stiff(2));
    M = M1 + M2;
    A = A1 + A2;
end

function [moment,area] = first_moment_of_area_rectangle(height,width,mean_height)
    area = height * width;
    moment = area * mean_height;
end

function [centroid, h_eq_over_h_3, height_max] = T_beam_properties(t_plate, h_stiff, width_plate, width_stiff)

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
    eta = eta_over_h * h; % location of neutral axis, from stiffened face of plate
    centroid = h - eta; % location of neutral axis, from unstiffened face of plate

    % eq 7.69 from MIT 2.080, but fixing typo: ^2 should be ^3 in term2,
    % and eta should be eta over h in term1
    % finding equivalent height that creates equal area moment of inertia
    term1 = 1 - 3 * eta_over_h + 3 * eta_over_h.^2;
    term2 = (H/h)^3 + 3*(H/h)^2 * eta_over_h + 3 * H/h * eta_over_h.^2;
    h_eq_over_h_3 = 4 * ( term1 + b./(2*a) .* term2 );

    height_max = h + H;

end
