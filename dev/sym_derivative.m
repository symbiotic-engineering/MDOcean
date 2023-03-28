syms a11 a12 a21 a22 b11 b12 b21 b22 x11 x12 x21 x22
A = [a11 a12; a21 a22];
B = [b11 b12; b21 b22];
X = [x11 x12; x21 x22];
eqn = A*X == B;
Xsol = solve(eqn,X)

n = 5;
syms Ann Bnn Xnn [n n]
eqn = Ann * Xnn == Bnn;
Xsol = solve(eqn,Xnn)
