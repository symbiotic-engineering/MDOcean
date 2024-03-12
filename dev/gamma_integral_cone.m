clear
syms sigma theta R a n d k
f_sigma = a * sigma + R;
C_cone_n(n) = (k * f_sigma * cos(theta))^n / factorial(n);
C_cone_sum = C_cone_n(0) - C_cone_n(2) + C_cone_n(4) - C_cone_n(6);
integrand = exp(k * sigma) * f_sigma * C_cone_sum;
integral = int( int(integrand, sigma, -d, 0), theta, 0, 2*pi);
integral = simplify(integral)