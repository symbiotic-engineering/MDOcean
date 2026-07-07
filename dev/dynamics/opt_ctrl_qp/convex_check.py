import cvxpy as cp
import numpy as np

def print_curvature(name, expr):
    print(f"\n=== {name} ===")
    print("Expression:", expr)
    print("is_affine   :", expr.is_affine())
    print("is_convex   :", expr.is_convex())
    print("is_concave  :", expr.is_concave())
    print("is_dcp      :", expr.is_dcp())
    print("is_quasiconvex :", expr.is_quasiconvex())
    print("is_quasiconcave:", expr.is_quasiconcave())


# ------------------------------------------------------------
# Variables (constraint limits)
# ------------------------------------------------------------
b_mu = cp.Variable(name="b_mu")
b_nu = cp.Variable(name="b_nu")

constraints = [b_mu >= 0, b_nu >= 0]

# ------------------------------------------------------------
# Closest-point candidate (one active constraint)
# Gamma_cl_sq = aff(b_mu) + sqrt(aff(b_mu))
# ------------------------------------------------------------
alpha = 2.0
beta = 1.0
aff_mu = alpha*b_mu + beta
constraints += [aff_mu >= 0]

Gamma_cl_sq = aff_mu + cp.sqrt(aff_mu)

print_curvature("Closest-point candidate: aff(b_mu) + sqrt(aff(b_mu))", Gamma_cl_sq)

# ------------------------------------------------------------
# Intersection candidate setup
#
# Construct Q(b_mu, b_nu) as a concave quadratic:
# Q = 5 - (b_mu - b_nu)^2
# ------------------------------------------------------------
Q = 5 - cp.square(b_mu - b_nu)
constraints += [Q >= 0]

print_curvature("Q = 5 - (b_mu - b_nu)^2", Q)
print_curvature("sqrt(Q)", cp.sqrt(Q))

# ------------------------------------------------------------
# Construct complex affine A(b_mu, b_nu)
#
# A = (A_r) + 1j*(A_i)
# where A_r and A_i are affine in b_mu, b_nu.
# ------------------------------------------------------------
A_r = 1.2*b_mu - 0.7*b_nu + 0.5
A_i = -0.4*b_mu + 0.9*b_nu + 0.1

A = A_r + 1j*A_i   # complex affine expression

print_curvature("Real part A_r", A_r)
print_curvature("Imag part A_i", A_i)
print_curvature("Complex affine A = A_r + 1j*A_i", A)

# ------------------------------------------------------------
# Form the complex intersection candidate:
# z = A + 1j*sqrt(Q)
# and compute |z|^2
#
# CVXPY represents |z|^2 as square(abs(z))
# ------------------------------------------------------------
z_plus = A + 1j*cp.sqrt(Q)
z_minus = A - 1j*cp.sqrt(Q)

Gamma_int_sq_plus = cp.square(cp.abs(z_plus))
Gamma_int_sq_minus = cp.square(cp.abs(z_minus))

print_curvature("|A + 1j*sqrt(Q)|^2", Gamma_int_sq_plus)
print_curvature("|A - 1j*sqrt(Q)|^2", Gamma_int_sq_minus)

# ------------------------------------------------------------
# Expand algebra manually to verify cross term:
#
# |A ± i*sqrt(Q)|^2 = A_r^2 + (A_i ± sqrt(Q))^2
#                   = A_r^2 + A_i^2 + Q ± 2*A_i*sqrt(Q)
# ------------------------------------------------------------
expanded_plus = cp.square(A_r) + cp.square(A_i) + Q + 2*A_i*cp.sqrt(Q)
expanded_minus = cp.square(A_r) + cp.square(A_i) + Q - 2*A_i*cp.sqrt(Q)

print_curvature("Expanded (+): A_r^2 + A_i^2 + Q + 2*A_i*sqrt(Q)", expanded_plus)
print_curvature("Expanded (-): A_r^2 + A_i^2 + Q - 2*A_i*sqrt(Q)", expanded_minus)

# ------------------------------------------------------------
# Compare the symbolic expressions (not numerically, but structurally)
# ------------------------------------------------------------
diff_plus = Gamma_int_sq_plus - expanded_plus
diff_minus = Gamma_int_sq_minus - expanded_minus

print_curvature("Difference: square(abs(z_plus)) - expanded_plus", diff_plus)
print_curvature("Difference: square(abs(z_minus)) - expanded_minus", diff_minus)

# ------------------------------------------------------------
# Try to build an optimization problem using these expressions
# ------------------------------------------------------------
objective = cp.Minimize(Gamma_int_sq_plus)
prob = cp.Problem(objective, constraints)

print("\n=== Problem check ===")
print("Problem is_dcp :", prob.is_dcp())
print("Problem is_dqcp:", prob.is_dqcp())