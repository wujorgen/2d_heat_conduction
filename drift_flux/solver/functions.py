import numpy as np
from constants import N_FACES, N_CELLS, dx, D, A_pipe, g_eff, dt
from helper_functions import fn_fric
from linear_solver import min_system


def calc_U_star_steady(p, rho, mu, u_last, u_in):
    def system(x):
        A = np.zeros((N_FACES, N_FACES))
        b = np.zeros((N_FACES))
        for i in range(N_FACES):
            if i == 0:  # inlet bc is velocity
                A[i, i] = 1
                b[i] = u_in
            elif i == N_FACES - 1:  # exit bc is pressure
                A[i, i] = 1
                A[i, i - 1] = - 1
                b[i] = 0
            else:
                A[i, i] = x[i] / dx
                A[i, i - 1] = - x[i] / dx
                rho_hat = rho[i:i+2].mean()  # avg i, i+1
                mu_hat = mu[i:i+2].mean()
                friction = fn_fric(rho_hat, x[i], mu_hat, D) / (2 * D) * x[i]**2
                b[i] = (- 1 / rho_hat) * (p[i+1] - p[i]) / dx - g_eff - friction
        return A, b
    sol = min_system(system, x0=u_last)
    return sol


def calc_U_star_transient(p, rho, mu, u_last):
    u_star = np.zeros(N_FACES)
    for i in range(1, N_FACES - 1):
        u_bar_ip1 = (u_last[i] + u_last[i+1]) / 2
        u_bar_i = (u_last[i] + u_last[i-1]) / 2
        convection_term = (u_bar_ip1**2 - u_bar_i**2) / dx
        rho_hat = rho[i:i+2].mean()  # avg i, i+1
        mu_hat = mu[i:i+2].mean()
        RHS = (-1 / rho_hat) * (p[i+1] - p[i]) / dx - g_eff - fn_fric(rho=rho_hat, v=u_last[i], mu=mu_hat, D=D) / (2 * D) * u_last[i]**2
        u_star[i] = u_last[i] + dt * (-convection_term + RHS)
    u_star[0] = u_star[1]
    u_star[-1] = u_star[-2]
    return {"x": u_star}


def calc_pressure_correction_steady(u, rho, dt=1):  # TODO: docstring. its called dt but its idfk what it is actually
    def system(x):
        A = np.zeros((N_CELLS, N_CELLS))
        b = np.zeros((N_CELLS))
        for i in range(N_CELLS):
            if i == 0:  # inlet bc is velocity -> no pressure gradient
                A[i, i] = 1
                A[i, i + 1] = -1
                b[i] = 0
            elif i == N_CELLS - 1:  # exit bc is pressure -> no correction 
                A[i, i] = 1
                b[i] = 0
            else:
                A[i, i + 1] = dt / dx**2
                A[i, i] = -2 * dt / dx**2
                A[i, i - 1] = dt / dx**2            

                rho_hat_iip1 = rho[i:i+2].mean()  # avg i, i+1
                rho_hat_iim1 = rho[i-1:i+1].mean()  # avg i-1, i
                b[i] = (rho_hat_iip1 * u[i] - rho_hat_iim1 * u[i-1]) / dx

        return A, b
    #print(*system(None))
    sol = min_system(system, x0=np.zeros(N_CELLS))
    return sol
    return {"x": np.linalg.solve(*system(None))}


def calc_pressure_correction_transient(u, rho, rho_last, dt=dt):
    def system(x):
        A = np.zeros((N_CELLS, N_CELLS))
        b = np.zeros((N_CELLS))
        for i in range(N_CELLS):
            if i == 0 or i == N_CELLS - 1:  # inlet bc is pressure -> no correction 
                A[i, i] = 1
                b[i] = 0
            elif i == 1 or i == N_CELLS - 2:  # inlet bc is pressure -> no correction 
                A[i, i] = 1
                b[i] = 0
            else:
                A[i, i + 1] = dt / dx**2
                A[i, i] = -2 * dt / dx**2
                A[i, i - 1] = dt / dx**2            

                rho_hat_iip1 = rho[i:i+2].mean()  # avg i, i+1
                rho_hat_iim1 = rho[i-1:i+1].mean()  # avg i-1, i
                b[i] = (rho_hat_iip1 * u[i] - rho_hat_iim1 * u[i-1]) / dx

                b[i] += (rho[i] - rho_last[i]) / dt

        return A, b
    #print(*system(None))
    sol = min_system(system, x0=np.zeros(N_CELLS))
    return sol


def gas_continuity_for_void_steady(rho_g, u_g, void_in, void_last, rho_g_last=None):
    def system(x):
        A = np.zeros((N_CELLS, N_CELLS))
        b = np.zeros((N_CELLS))
        for i in range(N_CELLS):
            if i == 0:
                # avg of ghost and first interior cell should be inlet void
                A[i, i] = 1/2
                A[i, i+1] = 1/2
                b[i] = void_in
            elif i == N_CELLS - 1:
                # no alpha gradient at exit
                A[i, i] = 1
                A[i, i-1] = -1
                b[i] = 0
            else:
                rho_hat_iip1 = rho_g[i:i+2].mean()  # avg i, i+1
                rho_hat_iim1 = rho_g[i-1:i+1].mean()  # avg i-1, i
                # A[i, i+1] = rho_hat_iip1 * u_g[i] / (2 * dx)
                # A[i, i  ] = rho_hat_iip1 * u_g[i] / (2 * dx) - rho_hat_iim1 * u_g[i-1] / (2 * dx)
                # A[i, i-1] = - rho_hat_iim1 * u_g[i-1] / (2 * dx)
                # upwind alpha based on velocity, instead of avg at faces (to avoid checkerboarding)
                idx_alpha_i = i if u_g[i] >= 0 else i + 1
                idx_alpha_im1 = i - 1 if u_g[i-1] >= 0 else i
                A[i, idx_alpha_i] += rho_hat_iip1 * u_g[i] / dx
                A[i, idx_alpha_im1] -= rho_hat_iim1 * u_g[i-1] / dx
                b[i] = 0
        return A, b
    sol = min_system(system, x0=void_last)
    return sol


def gas_continuity_for_void_transient(rho_g, u_g, void_in, void_last, rho_g_last):
    """
        The equations are discretized from cell face to cell face.

    :param void_last: alpha at last time step
    :param rho_g_last: gas density at last time step
    """
    def system(x):
        A = np.zeros((N_CELLS, N_CELLS))
        b = np.zeros((N_CELLS))
        for i in range(N_CELLS):
            if i == 0 or i == 1:
                # avg of ghost and first interior cell should be inlet void
                #A[i, i] = 1/2
                #A[i, i+1] = 1/2
                A[i, i] = 1
                b[i] = void_in
            elif i == N_CELLS - 1:
                # no alpha gradient at exit
                A[i, i] = 1
                A[i, i-1] = -1
                b[i] = 0
            else:
                rho_hat_iip1 = rho_g[i:i+2].mean()  # avg i, i+1
                rho_hat_iim1 = rho_g[i-1:i+1].mean()  # avg i-1, i
                # A[i, i+1] = rho_hat_iip1 * u_g[i] / (2 * dx)
                # A[i, i  ] = rho_hat_iip1 * u_g[i] / (2 * dx) - rho_hat_iim1 * u_g[i-1] / (2 * dx)
                # A[i, i-1] = - rho_hat_iim1 * u_g[i-1] / (2 * dx)
                # upwind alpha based on velocity, instead of avg at faces (to avoid checkerboarding)
                idx_alpha_i = i if u_g[i] >= 0 else i + 1
                idx_alpha_im1 = i - 1 if u_g[i-1] >= 0 else i
                A[i, idx_alpha_i] += dt * rho_hat_iip1 * u_g[i] / dx
                A[i, idx_alpha_im1] -= dt * rho_hat_iim1 * u_g[i-1] / dx
                A[i, i] += rho_g[i]
                b[i] = void_last[i] * rho_g_last[i]
        return A, b
    sol = min_system(system, x0=void_last)
    return sol