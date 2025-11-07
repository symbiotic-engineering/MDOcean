%% symbolically
b_0_symbolic

%% functions
function b_0_symbolic()
    % from A_b_c_matrix_N10_M10_K10_heaving_outer, tracing back the top (k=0)
    % term of the bottom block of the b-vector
    clear
    syms a1 a2 d1 d2 h m0 real positive

    t16 = d2.*2.0;
    t19 = h.*2.0;
    t35 = m0.*t19;
    t51 = -h;
    t52 = -t19;
    t53 = 1.0./h;
    t54 = 1.0./m0;
    t81 = sinh(t35);
    t161 = d2+t51;
    t244 = m0.*t161;
    t255 = t16+t52;
    t261 = sinh(t244);
    t295 = 1.0./t255;

    b0 = -a2.*t54.*t261.*t295.*1.0./sqrt((t53.*t54.*t81)./4.0+1.0./2.0);
    b0 = subs(expand(b0),sinh(m0*h), cosh(m0*h)*tanh(m0*h));
    [num,den] = numden(simplifyFraction(b0));
    new_num = simplifyFraction(num/cosh(h*m0));
    new_den = sqrt(partfrac(den^2/cosh(h*m0)^2));
    pretty(new_num/new_den)
    b0_approx = subs(new_num/new_den, {tanh(h*m0),1/cosh(h*m0)^2},{1,0});
    pretty(simplify(b0_approx))
end
