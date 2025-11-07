clear;
syms sigma theta R a n d k;
f_sigma = a * sigma + R;
C_cone_n(n) = (k * f_sigma * cos(theta))^(2 * n) / factorial(2 * n) * (-1)^n;
C_cone_sum = symsum(C_cone_n, n, 0, Inf);
integrand = exp(k * sigma) * f_sigma * C_cone_sum;
integral = int(int(integrand, sigma, -d, 0), theta, 0, 2 * pi);
integral = simplify(integral);
