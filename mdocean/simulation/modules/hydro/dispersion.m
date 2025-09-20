function k = dispersion(w,h,g)
    h_lambda_deep = 0.4; % h/lambda > threshold for deep
    h_lambda_shallow = 0.05; % h/lambda < threshold for shallow
    idx_deep = w > sqrt(h_lambda_deep * 2*pi * g / h); 
    idx_shallow = w < h_lambda_shallow * 2*pi * sqrt(g/h);
    idx_mid = ~idx_deep & ~idx_shallow;

    k = zeros(size(w));
    k(idx_deep) = w(idx_deep).^2 / g;
    k(idx_shallow) = w(idx_shallow) / sqrt(g*h);

    if any(idx_mid,'all')
        err = 1;
        iters = 0;
        k_deep_guess = w(idx_mid).^2 / g;
        k_guess = k_deep_guess;
        while err > 0.2 && iters < 50
            % fixed point iteration, dividing by k_deep_guess for better stability
            k_new = w(idx_mid).^2 / g ./ tanh(k_guess*h);
            err = max(abs(k_new - k_guess)./k_deep_guess,[],'all');
            k_guess = k_new;
            iters = iters + 1;
        end
        k(idx_mid) = k_new;
    end
end