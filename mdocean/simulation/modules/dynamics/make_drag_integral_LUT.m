function [drag_integral_fcn] = make_drag_integral_LUT(rp_min, rp_step, rp_max, theta_step, kappa_step, kappa_max)
    % r_prime: magnitude of wave velocity over magnitude of wec velocity (r in paper)
    r_prime_vec = [rp_min : rp_step : rp_max 1e6];
    r_vec = 1./r_prime_vec;
    theta_vec = -pi/2 : theta_step : pi/2;
    kappa_vec = 0 : kappa_step : kappa_max;

    [R,TH,K] = meshgrid(r_vec,theta_vec,kappa_vec);

    [B,G_real,G_imag] = compute_drag_integral(R,TH,K);

    drag_integral_LUT_fcn = ...
        @(r_q,th_q,k_q) deal(interp3(R,TH,K,B,     r_q,wrapToPiOver2(th_q),k_q),...
                             interp3(R,TH,K,G_real,r_q,wrapToPiOver2(th_q),k_q),...
                             interp3(R,TH,K,G_imag,r_q,wrapToPiOver2(th_q),k_q));

    drag_integral_fcn = @(r,th,k) drag_LUT_or_analytical_or_compute(r,th,k,drag_integral_LUT_fcn);
end

function [B,G_real,G_imag] = drag_LUT_or_analytical_or_compute(R,TH,K,LUT_fcn)
    [B,G_real,G_imag] = deal(zeros(size(R)));
    
    % first compute analytically where possible (only r_prime=0, R=0, and K=0 limit cases)
    rp = 1./R;
    idx_rp_zero = rp < 1e-6;
    B(idx_rp_zero) = pi/2;
    G_imag(idx_rp_zero) = 0;
    G_real(idx_rp_zero) = pi * besselj(1,K(idx_rp_zero)) ./ (K(idx_rp_zero));

    idx_r_zero = R == 0 & K ~= 0;
    B(idx_r_zero) = 0;
    G_imag(idx_r_zero) = 0;
    G_real(idx_r_zero) = pi * besselj(1,K(idx_r_zero)) ./ (K(idx_r_zero));

    idx_k_zero = K == 0 & R ~= 0;
    k_zero_const_term = pi/2 * sqrt(1 + rp(idx_k_zero).^2 - 2*rp(idx_k_zero).*sin(TH(idx_k_zero)));
    B(idx_k_zero) = k_zero_const_term;
    G_imag(idx_k_zero) = 0;
    G_real(idx_k_zero) = k_zero_const_term;

    idx_r_k_zero = R == 0 & K == 0;
    B(idx_r_k_zero) = 0;
    G_imag(idx_r_k_zero) = 0;
    G_real(idx_r_k_zero) = pi/2;

    % then use LUT everywhere else (queries outside the LUT will return NaN)
    idx_LUT = ~idx_rp_zero & ~idx_k_zero & ~idx_r_zero & ~idx_r_k_zero;
    [B(idx_LUT), G_real(idx_LUT), G_imag(idx_LUT)] = LUT_fcn(R(idx_LUT),TH(idx_LUT),K(idx_LUT));

    % for NaNs outside the LUT, compute directly (this is the slowest but should be rare)
    idx_nan = isnan(B) & ~isnan(R);
    if any(idx_nan)
         disp(['Computing drag integral numerically for rp=' mat2str(rp(idx_nan)) ', kappa=' mat2str(K(idx_nan)) ...
         ' because they are outside the range of the precomputed lookup table. If you see this message often' ...
         ' and want to reduce compute time, go to parameters.m, find the call to make_drag_integral_LUT(),'... 
         ' and adjust the rp_min, rp_max, and/or kappa_max inputs to expand the lookup table.'])
         [B(idx_nan), G_real(idx_nan), G_imag(idx_nan)] = compute_drag_integral(R(idx_nan),TH(idx_nan),K(idx_nan));
    end
    
end

function out = wrapToPiOver2(angle)

    wrapped = wrapToPi(angle);
    out = sign(wrapped) .* min( abs(wrapped), pi - abs(wrapped));
end
