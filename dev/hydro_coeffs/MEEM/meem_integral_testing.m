syms r r1 r2 real positive;
syms a real;
integrand = @(r) r.*besselk(0,r*a);
int(integrand,r,r1,r2)