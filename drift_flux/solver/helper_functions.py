import numpy as np


# optimization routine for implicit solves
# (quadratic) obj: min ||Ax-b||_2^2
fn_f = lambda A, x, b: x.T @ A.T @ A @ x - 2 * b.T @ A @ x + b @ b

# gradient and hessian
fn_g = lambda A, x, b: 2 * A.T @ (A @ x - b)
fn_h = lambda A, x, b: 2 * A.T @ A


def fn_fric(rho, v, mu, D):
    """Calculates Blasius friction factor."""
    Re_m = rho * v * D / mu
    fB = 0.316 / (np.abs(Re_m+0.001)**0.25)
    # print(Re_m, rho, v, D, mu)
    return fB * np.sign(Re_m)


def fn_check_step_convergence(sol:dict, itr, name):
    if "converged" in sol.keys():
        if not sol["converged"]:
            print(f"Warning at iteration {itr}: {name} not converged.")