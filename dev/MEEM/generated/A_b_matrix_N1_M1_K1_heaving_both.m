function [A,b] = A_b_matrix_N1_M1_K1_heaving_both(a1,a2,d1,d2,h,m0,m_k_const1)
%A_b_matrix_N1_M1_K1_heaving_both
%    [A,B] = A_b_matrix_N1_M1_K1_heaving_both(A1,A2,D1,D2,H,M0,M_K_CONST1)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    19-May-2024 17:28:10

t2 = a2.*m0;
t3 = a2.*m_k_const1;
t4 = pi.^2;
t5 = d1.*2.0;
t6 = d2.*2.0;
t7 = d1.^2;
t8 = d2.^2;
t9 = h.*2.0;
t10 = h.^2;
t11 = m0.^2;
t12 = m_k_const1.^2;
t17 = 1.0./pi;
t19 = 1.0./a2;
t20 = -d2;
t21 = -h;
t23 = 1.0./h;
t24 = 1.0./m0;
t25 = 1.0./m_k_const1;
t26 = d2.*h.*-2.0;
t27 = sqrt(2.0);
t30 = d1./2.0;
t31 = h./2.0;
t32 = sqrt(h);
t33 = m0.^(3.0./2.0);
t38 = 1.0./sqrt(m0);
t13 = h.*t5;
t14 = h.*t6;
t15 = m0.*t9;
t16 = m_k_const1.*t9;
t18 = 1.0./t4;
t22 = -t9;
t34 = -t7;
t35 = d1+t21;
t36 = d2+t21;
t37 = h+t20;
t39 = t8.*t11;
t40 = t8.*t12;
t41 = t10.*t11;
t42 = t10.*t12;
t52 = t11.*t26;
t28 = sinh(t15);
t29 = sin(t16);
t43 = t11.*t14;
t44 = t12.*t14;
t45 = m0.*t36;
t46 = m_k_const1.*t36;
t47 = t5+t22;
t48 = t6+t22;
t49 = t35.^2;
t50 = t35.^3;
t51 = t36.^2;
t55 = 1.0./t35;
t56 = 1.0./t36;
t57 = -t40;
t58 = -t42;
t81 = t8+t13+t26+t34;
t84 = -1.0./(t7-t8-t13+t14);
t86 = t4+t39+t41+t52;
t53 = sinh(t45);
t54 = sin(t46);
t59 = t15+t28;
t60 = 1.0./t47;
t61 = 1.0./t48;
t62 = a1.*t55.*pi;
t63 = a1.*t56.*pi;
t64 = a2.*t55.*pi;
t65 = a2.*t56.*pi;
t76 = t35.*t56.*pi;
t77 = (t23.*t25.*t29)./2.0;
t87 = 1.0./t86;
t88 = t4+t44+t57+t58;
t66 = besseli(0,t63);
t67 = besseli(0,t64);
t68 = besseli(0,t65);
t69 = besseli(1,t65);
t70 = -t63;
t71 = -t65;
t72 = 1.0./sqrt(t59);
t80 = sin(t76);
t83 = t77+1.0;
t89 = 1.0./t88;
t73 = besselk(0,t70);
t74 = besselk(0,t71);
t75 = besselk(1,t71);
t78 = 1.0./t67;
t79 = 1.0./t68;
t85 = 1.0./sqrt(t83);
t82 = 1.0./t74;
mt1 = [-t30+t31,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t35.*t78.*besseli(0,t62),0.0,0.0,0.0,(t35.*t36.*t78.*t80.*besseli(1,t62).*-2.0)./(t7-t8-t13+t14),0.0,0.0,t30-t31,0.0,d2.*(-1.0./2.0)+t31,0.0,0.0,0.0,0.0,0.0,t17.*t27.*t36.*t66.*t79.*t80,(t17.*t36.*t49.*t66.*t79.*t80.*-2.0)./(t7-t8-t13+t14),0.0,t37,0.0,-t79.*pi.*besseli(1,t63),t27.*t32.*t33.*t36.*t53.*t69.*t72.*t79.*t87.*pi.*-2.0,t46.*t54.*t69.*t79.*t85.*t89.*pi.*2.0,(t35.*log(a1.*t19))./2.0,0.0,0.0,0.0,(t36.*(-1.0./2.0))./a1,0.0,t19.*t32.*t38.*t53.*t72,(t27.*t54.*t85)./(t3.*2.0),t17.*t27.*t36.*t73.*t80.*t82];
mt2 = [(t17.*t36.*t49.*t73.*t80.*t82.*-2.0)./(t7-t8-t13+t14),0.0,t37,0.0,-t82.*pi.*besselk(1,t70),t27.*t32.*t33.*t36.*t53.*t72.*t75.*t82.*t87.*pi.*-2.0,t46.*t54.*t75.*t82.*t85.*t89.*pi.*2.0,0.0,0.0,t32.*t38.*t53.*t72.*2.0,t27.*t32.*t33.*t51.*t53.*t72.*t87.*-2.0,0.0,0.0,(m0.*t21.*besselh(1.0,1,t2))./besselh(0.0,1,t2),0.0,0.0,0.0,t25.*t27.*t54.*t85,t36.*t46.*t54.*t85.*t89.*2.0,0.0,0.0,0.0,(m_k_const1.*t21.*besselk(1,t3))./besselk(0,t3)];
A = reshape([mt1,mt2],8,8);
if nargout > 1
    b = [(t56.*(d1+t20).*(t7.*2.0-d1.*h.*4.0+h.*t9-a1.^2.*3.0))./1.2e+1;t18.*t27.*t50.*t60.*2.0-t18.*t27.*t50.*t61.*2.0;t8.*(-1.0./6.0)-t10./6.0+(d2.*h)./3.0+1.0./t19.^2./4.0;t18.*t27.*1.0./t56.^3.*t61.*2.0;0.0;-a1.*t17.*t27.*t36.*t60.*t80;-a2.*t24.*t53.*t61.*1.0./sqrt((t23.*t24.*t28)./4.0+1.0./2.0);-a2.*t25.*t54.*t61.*1.0./sqrt((t23.*t25.*t29)./4.0+1.0./2.0)];
end