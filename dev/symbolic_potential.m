syms r theta z real
syms a_n k0 K lambda_n h g A omega m n b1 b2 h1 h2 real positive
syms C2_mn_B C2_mn_C C_n_D C_n_B C_n_C C_m_B C_m_C C_mn_D C_mn_B C_mn_C

epsilon_m = 2;
C_m_I = epsilon_m * 1i^m;

R_mn_D = besselk(m,a_n*r);
R_mn_B = besselk(m,a_n*r) + C2_mn_B * besseli(m,a_n*r);
R_mn_C = besselk(m,a_n*r) + C2_mn_C * besseli(m,a_n*r);
R_m_I  = besselj(m,k0* r);

Z_n_D = C_n_D * cosh(k0 * z) + cos(lambda_n*z);
Z_n_B = C_n_B * cosh(k0 * z) + cos(lambda_n*z);
Z_n_C = C_n_C * cosh(k0 * z) + cos(lambda_n*z);
Z_I = cosh(k0*z) / cosh(k0 * h);

theta_m_D = exp(1i*m*theta);
theta_m_B = exp(1i*m*theta) + C_m_B * exp(-1i*m*theta);
theta_m_C = exp(1i*m*theta) + C_m_C * exp(-1i*m*theta);
theta_m_I = cos(m*theta);

phiD_mn = C_mn_D * R_mn_D * Z_n_D * theta_m_D;
phiB_mn = C_mn_B * R_mn_B * Z_n_B * theta_m_B;
phiC_mn = C_mn_C * R_mn_C * Z_n_C * theta_m_C;
phiI_m  = C_m_I  * R_m_I  * Z_I   * theta_m_I;

const = - 1i * g * A / omega;
phiD = const * symsum(symsum(phiD_mn,n,0,Inf),m,-Inf,Inf);
phiB = const * symsum(symsum(phiB_mn,n,0,Inf),m,-Inf,Inf);
phiC = const * symsum(symsum(phiC_mn,n,0,Inf),m,-Inf,Inf);
phiI = const * (symsum(phiI_m,m,-Inf,-1) + 1/2*symsum(phiI_m,m,0,0) + symsum(phiI_m,m,1,Inf));
phiA = phiD + phiI;

dphiA = jacobian(phiA,[r z]);
dphiB = jacobian(phiB,[r z]);
dphiC = jacobian(phiC,[r z]);

BC1 = -K * phiA + subs(dphiA(2),z,h) == 0;
BC2 = subs(dphiA(1),r,b2) == 0;
BC3 = subs(dphiB(1),r,b1) == 0;
BC4 = subs(dphiA(2),z,h2) == 0;
BC5 = subs(dphiC(2),z,h1) == 0;
BC6 = subs(phiA,r,b2) == subs(phiB,r,b2);
BC7 = subs(phiC,r,b1) == subs(phiB,r,b1);
BC8 = subs(dphiA(1),r,b2) == subs(dphiB(1),r,b2);
BC9 = subs(dphiC(1),r,b1) == subs(dphiB(1),r,b1);
BC11 = subs(dphiA(2),z,0) == 0;
BC12 = subs(dphiB(2),z,0) == 0;
BC13 = subs(dphiC(2),z,0) == 0;

eqns = [BC1 BC2 BC3 BC4 BC5 BC6 BC7 BC8 BC9 BC11 BC12 BC13];
vars = [C2_mn_B C2_mn_C C_n_D C_n_B C_n_C C_m_B C_m_C C_mn_D C_mn_B C_mn_C];
%solve(eqns, vars)

