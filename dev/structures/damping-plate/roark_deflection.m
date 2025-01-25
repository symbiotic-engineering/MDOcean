syms r0_point r0_dist w a b v q E t

C2 = 1/4*(1-(b/a)^2*(1+2*log(a/b)));
C3 = b/4/a*(((b/a)^2+1)*log(a/b)+(b/a)^2-1);
C8 = 1/2*(1+v+(1-v)*(b/a)^2);
C9 = b/a*((1+v)/2*log(a/b)+(1-v)/4*(1-(b/a)^2));
L3 = r0_point/4/a*(((r0_point/a)^2+1)*log(a/r0_point)+(r0_point/a)^2-1);
L9 = r0_point/a*((1+v)/2*log(a/r0_point)+(1-v)/4*(1-(r0_point/a)^2));
L11 = 1/64*(1+4*(r0_dist/a)^2-5*(r0_dist/a)^4-4*(r0_dist/a)^2*(2+(r0_dist/a)^2)*log(a/r0_dist));
L17 = 1/4*(1-(1-v)/4*(1-(r0_dist/a)^4)-(r0_dist/a)^2*(1+(1+v)*log(a/r0_dist)));

Mrb = -q*a^2/C8 * (C9*(a^2-r0_dist^2)/(2*a*b) - L17);
Qb = q/2/b * (a^2 - r0_dist^2);
D = E*t^3/12/(1-v^2);

l1 = -w*a^3/D * (C2/C8 * (r0_point*C9/b - L9) - (r0_point*C3/b) + L3);
l2 = Mrb*a^2*C2/D + Qb*a^3*C3/D - q*a^4*L11/D;
y = l1 + l2;

matlabFunction(y,l1,l2, 'File','Roark_func','Vars', [r0_point,r0_dist,a,w,t,E,q,b,v]);