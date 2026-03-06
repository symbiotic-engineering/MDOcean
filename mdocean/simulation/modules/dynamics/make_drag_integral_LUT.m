function [drag_integral_LUT_fcn] = make_drag_integral_LUT(rp_min, rp_step, rp_max, theta_step, kappa_step, kappa_max)
    r_primve_vec = rp_min : rp_step : rp_max;
    r_vec = 1./r_primve_vec;
    theta_vec = -pi/2 : theta_step : pi/2;
    kappa_vec = 0 : kappa_step : kappa_max;

    [R,TH,K] = meshgrid(r_vec,theta_vec,kappa_vec);

    T = compute_drag_integral(R,TH,K);

    drag_integral_LUT_fcn = @(r_q,th_q,k_q) interp3(R,TH,K,T,r_q,th_q,k_q);
end

