% from A_b_c_matrix_N10_M10_K10_heaving_outer, tracing back the top (k=0)
% term of the bottom block of the b-vector
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
pretty(b0)