function [V_g, mult] = group_velocity(w, k, g, h)

    % notebook p129 2/25/25
    % w2_over_gk and mult equal 1 in deep water, not 1 in finite depth
    w2_over_gk = w.^2 ./ (g * k);
    mult = k * h .* (1 - w2_over_gk.^2) + w2_over_gk;

    % group velocity: defined as dw/dk
    V_g = g ./ (2 * w) .* mult;

end
