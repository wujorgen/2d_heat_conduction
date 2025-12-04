import numpy as np
from helper_functions import fn_f, fn_g, fn_h

def SPAI(A, m):
    """Perform m step of the SPAI iteration."""
    from scipy.sparse import identity
    from scipy.sparse import diags
    from scipy.sparse.linalg import onenormest
    
    n = A.shape[0]
    
    ident = identity(n, format='csr')
    alpha = 2 / onenormest(A @ A.T)
    M = alpha * A
        
    for index in range(m):
        C = A @ M
        G = ident - C
        AG = A @ G
        trace = (G.T @ AG).diagonal().sum()
        alpha = trace / np.linalg.norm(AG.data)**2
        M = M + alpha * G
        
    return np.array(M)


def min_system(system, x0, lr = 0.5, tol_itr = 1e-7, tol_grad = 1e-4, max_itr:int = int(1e3), verbose:bool = False) -> dict:
    """Solves a linear system (Ax=b) by minimizing a quadratic objective function ||Ax-b||_2^2

    :param callable system: provides the coefficient matrix and RHS of the system.
                    must be provided as function, assuming that A and b are functions of x.
    :param int size: size of the solution vector x
    :return sol: dictionary 
    """
    f_log = []
    g_log = []
    h_log = []
    x_log = []
    r_log = []

    #x0 = np.full(NX, 1)  # TODO: can't start this at zero, or too small of a value.
    x = x0
    converged = False

    for itr in range(max_itr):
        A, b = system(x)
        # M = SPAI(A, 50)  # for preconditioning
        # A = M @ A
        # b = M @ b

        f = fn_f(A, x, b); f_log.append(f)
        g = fn_g(A, x, b); g_log.append(g)
        h = fn_h(A, x, b); h_log.append(h)

        try:
            step = lr * np.linalg.inv(h) @ g  # newton step (for minimizing)
        except:
            breakpoint()
        # step = lr * g  # gradient step

        x_new = x - step

        r = np.sqrt(np.sum(np.mean((x_new - x)**2))); r_log.append(r)
        gradnorm = np.linalg.norm(g)

        if r <= tol_itr:
            if verbose:
                print("Solution converged by residual tolerance at iteration", itr)
            converged = True
            break
        if gradnorm <= tol_grad:
            if verbose:
                print("Solution converged by gradient tolerance at iteration", itr)
            converged = True
            break

        x = x_new  # if not converged, prepare for next iteration
    
    return {
        "x": x,
        "converged": converged,
        "f_log": f_log,
        "g_log": g_log,
        "r_log": r_log,
    }
