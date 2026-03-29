function [drag_integral_fcn] = make_drag_integral_LUT(rp_min, rp_step, rp_max, theta_step, kappa_step, kappa_max)
    % r_prime: magnitude of wec velocity over magnitude of wave velocity (1/r in paper), used for polar plot radius
    r_prime_vec = [1e-6 rp_min : rp_step : rp_max 1e6];
    r_vec = 1./r_prime_vec;
    theta_vec = -pi/2 : theta_step : pi/2;
    kappa_vec = 0 : kappa_step : kappa_max;

    [R,TH,K] = meshgrid(fliplr(r_vec),theta_vec,kappa_vec);

    [B,G_real,G_imag] = compute_drag_integral(R,TH,K);

    F_B      = myinterpFcn(R,TH,K,B);
    F_G_real = myinterpFcn(R,TH,K,G_real);
    F_G_imag = myinterpFcn(R,TH,K,G_imag);

    drag_integral_LUT_fcn = @(r_q,th_q,k_q) permuteInterpFcn(r_q,th_q,k_q,F_B,F_G_real,F_G_imag);

%     p = [2  1 3];
%     drag_integral_LUT_fcn_2 = ...
%         @(r_q,th_q,k_q) deal( permute(      F_B( permute(r_q,p), wrapToPiOver2(permute(th_q,p)), permute(k_q,p) ),p),...
%                               permute( F_G_real( permute(r_q,p), wrapToPiOver2(permute(th_q,p)), permute(k_q,p) ),p),...
%                               permute( F_G_imag( permute(r_q,p), wrapToPiOver2(permute(th_q,p)), permute(k_q,p) ),p) );
%     test_in = [1 2; 3 4]/4;
%     test_in(:,:,2) = [5 6; 7 8]/4;
%     [a, b, c] = drag_integral_LUT_fcn(  test_in,2*test_in,3*test_in);
%     [a2,b2,c2]= drag_integral_LUT_fcn_2(test_in,2*test_in,3*test_in);

    drag_integral_fcn = @(r,th,k) drag_LUT_or_analytical_or_compute(r,th,k,drag_integral_LUT_fcn);
end

function [out_1,out_2,out_3] = permuteInterpFcn(r_q,th_q,k_q,F1,F2,F3)
    p = [2 1 3];
    r_p = permute(r_q,p);
    th_p = wrapToPiOver2(permute(th_q,p));
    k_p = permute(k_q,p);

    out_1_p = F1(r_p,th_p,k_p);
    out_2_p = F2(r_p,th_p,k_p);
    out_3_p = F3(r_p,th_p,k_p);

    out_1 = permute(out_1_p,p);
    out_2 = permute(out_2_p,p);
    out_3 = permute(out_3_p,p);

end

function [B,G_real,G_imag] = drag_LUT_or_analytical_or_compute(R,TH,K,LUT_fcn)
    B      = zeros(size(R));
    G_real = zeros(size(R));
    G_imag = zeros(size(R));
    
    R(R > 1e6) = Inf;
    R(R < 1e-6) = 0;

    % first compute analytically where possible (only R=Inf, R=0, and K=0 limit cases)

    % r -> Inf
    idx_r_inf = isinf(R) & K ~= 0;
    B(idx_r_inf) = 0;
    G_r_inf = pi * besselj(1,K(idx_r_inf)*(1i-1)) ./ (K(idx_r_inf)*(1i-1));
    G_imag(idx_r_inf) = imag(G_r_inf);
    G_real(idx_r_inf) = real(G_r_inf);

    idx_r_inf_k_zero = isinf(R) & K == 0;
    B(idx_r_inf_k_zero) = 0;
    G_imag(idx_r_inf_k_zero) = 0;
    G_real(idx_r_inf_k_zero) = pi/2;

    % r -> 0
    idx_r_zero = R == 0 & K ~= 0;
    B(idx_r_zero) = pi/2;
    G_imag(idx_r_zero) = 0;
    G_real(idx_r_zero) = pi * besselj(1,K(idx_r_zero)) ./ (K(idx_r_zero));

    idx_r_zero_k_zero = R == 0 & K == 0;
    B(idx_r_zero_k_zero) = pi/2;
    G_imag(idx_r_zero_k_zero) = 0;
    G_real(idx_r_zero_k_zero) = pi/2;

    % k -> 0
    idx_k_zero = K == 0 & ~isinf(R) & R ~= 0;
    k_zero_const_term = pi/2 * sqrt(1 + R(idx_k_zero).^2 - 2*R(idx_k_zero).*sin(TH(idx_k_zero)));
    B(idx_k_zero) = k_zero_const_term;
    G_imag(idx_k_zero) = 0;
    G_real(idx_k_zero) = k_zero_const_term;


    % then use LUT everywhere else (queries outside the LUT will return NaN)
    idx_LUT = ~idx_r_inf & ~idx_r_inf_k_zero & ~idx_r_zero & ~idx_r_zero_k_zero & ~idx_k_zero;
    [B(idx_LUT), G_real(idx_LUT), G_imag(idx_LUT)] = LUT_fcn(R(idx_LUT),TH(idx_LUT),K(idx_LUT));

    % for NaNs outside the LUT, compute directly (this is the slowest but should be rare)
    idx_nan = isnan(B) & ~isnan(R);
    if any(idx_nan)
         disp(['Computing drag integral numerically for rp=' mat2str(1./R(idx_nan)) ', kappa=' mat2str(K(idx_nan)) ...
         ' because they are outside the range of the precomputed lookup table. If you see this message often' ...
         ' and want to reduce compute time, go to calkit.yaml, find the call to make_drag_integral_LUT(),'... 
         ' and adjust the rp_min, rp_max, and/or kappa_max inputs to expand the lookup table.'])
         [B(idx_nan), G_real(idx_nan), G_imag(idx_nan)] = compute_drag_integral(R(idx_nan),TH(idx_nan),K(idx_nan));
    end
    
end

function out = wrapToPiOver2(angle)

    wrapped = wrapToPi(angle);
    out = sign(wrapped) .* min( abs(wrapped), pi - abs(wrapped));
end

function F = myinterpFcn(X,Y,Z,V)
    p = [2 1 3];
    X = permute(X,p);
    Y = permute(Y,p);
    Z = permute(Z,p);
    V = permute(V,p);
    F = griddedInterpolant(X,Y,Z,V,'linear','none');
end
