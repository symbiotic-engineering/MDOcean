function [X_slam_min, X_slam_max,...
          idx_imag, diameter_margin, theta] = get_slamming_min_max(A, k, D, phase_X, delta_z)

    idx_large_wave = A > delta_z;

    theta_small_wave = pi + max(0, -k * D / 2 + abs(pi - phase_X));
    theta_large_wave = phase_X + k * D / 2 .* sign(pi - phase_X);
    theta = theta_small_wave;
    theta(idx_large_wave) = theta_large_wave(idx_large_wave);

    A_term = A .* cos(theta);
    sqrt_term = sqrt( delta_z^2 + A_term.^2 - A.^2 );
    
    X_slam_max =  sqrt_term + A_term;
    X_slam_min = -sqrt_term + A_term;

    idx_imag = imag(sqrt_term)~=0; % case where slamming occurs even for stationary body
    
    if any(idx_large_wave(:))
        k_max_large_wave = max(k(idx_large_wave),[],'all');
        A_max_large_wave = max(A(idx_large_wave),[],'all');
        diameter_margin = asin(delta_z/A_max_large_wave) - k_max_large_wave * D/2;
    else
        diameter_margin = 1; % always ok if no large waves
    end

end