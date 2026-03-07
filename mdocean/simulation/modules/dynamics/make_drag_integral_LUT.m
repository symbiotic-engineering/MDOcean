function [drag_integral_LUT_fcn] = make_drag_integral_LUT(rp_min, rp_step, rp_max, theta_step, kappa_step, kappa_max)
    r_primve_vec = [rp_min : rp_step : rp_max 1e6];
    r_vec = 1./r_primve_vec;
    theta_vec = -pi/2 : theta_step : pi/2;
    kappa_vec = 0 : kappa_step : kappa_max;

    [R,TH,K] = meshgrid(r_vec,theta_vec,kappa_vec);

    [B,G_real,G_imag] = compute_drag_integral(R,TH,K);

    drag_integral_LUT_fcn = ...
        @(r_q,th_q,k_q) deal(interp3(R,TH,K,B,     r_q,wrapToPiOver2(th_q),k_q),...
                             interp3(R,TH,K,G_real,r_q,wrapToPiOver2(th_q),k_q),...
                             interp3(R,TH,K,G_imag,r_q,wrapToPiOver2(th_q),k_q));
end

function out = wrapToPiOver2(angle)

    wrapped = wrapToPi(angle);
    out = sign(wrapped) .* min( abs(wrapped), pi - abs(wrapped));
end
