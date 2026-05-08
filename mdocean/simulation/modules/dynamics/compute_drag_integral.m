function [B,G_real,G_imag] = compute_drag_integral(R,TH,K)
    [B,G_real,G_imag] = deal(zeros(size(R)));

    for i=1:numel(R)
        B(i)      = integral(@(x) f(x, R(i), TH(i), K(i)),               -1, 1);
        G_real(i) = integral(@(x) f(x, R(i), TH(i), K(i)).*cos(-K(i)*x), -1, 1);
        G_imag(i) = integral(@(x) f(x, R(i), TH(i), K(i)).*sin(-K(i)*x), -1, 1);
    end
end

function f = f(x,r,theta,kappa)
    t1 = 1 - r .* sin(theta) .* exp(-kappa .* x);
    t2 =     r .* cos(theta) .* exp(-kappa .* x);
    f = sqrt(1-x.^2) .* sqrt( t1.^2 + t2.^2 );
end