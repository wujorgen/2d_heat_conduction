# %%
import numpy as np
import scipy
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from dataclasses import dataclass
from copy import deepcopy

# Finite Difference Grid
NX = 201
L = 2
X = np.linspace(0, L, NX, endpoint=True)
dx = L / (NX - 1)
D = 1 / 39.37
dt = 0.01
T_end = 1
NT = T_end / dt + 1

# orientation
theta = 90  # degrees
g_constant = 9.81
g_eff = g_constant * np.sin(np.deg2rad(theta))

# Drift flux parameters
C0 = 1.2
V_gj = 0.35 * np.sqrt(g_eff * D)

# Properties for air/water system

T_AMB = 273.15 + 20  # room temp?
P_ATM = 101325  # N/m^2

air = 'air'

# 1. Get the universal gas constant (J/mol/K)
R = PropsSI('gas_constant', air)

# 2. Get the molar mass (kg/mol)
M = PropsSI('molemass', air) # in kg/mol (Note: CoolProp uses kg/mol by default for this)

# 3. Calculate the specific gas constant (J/kg/K)
Rbar = R / M

fn_rho_g = lambda P, T: P / (Rbar * T)
mu_g = PropsSI("V", "P", P_ATM, "T", T_AMB, "Air")
rho_l = PropsSI("D", "P", P_ATM, "T", T_AMB, "Water")
mu_l = PropsSI("V", "P", P_ATM, "T", T_AMB, "Water")

print("Air Properties")
print(f"Universal gas constant (R): {R:.3f} J/(mol·K)")
print(f"Molar mass (M): {M:.3f} kg/mol")
print(f"Specific gas constant (Rs): {Rbar:.3f} J/(kg·K)")
print(f"Air density at room temp should be around 1.204 kg/m^3. Confirm: {fn_rho_g(P_ATM, T_AMB):.3f} kg/m^3")
print(f"Water density at room temp: {rho_l:.3f} kg/m^3")
print(f"Water dynamic viscosity at room temp: {mu_l:.3f} kg/m^3")


# inlet velocity boundary condition
INLET_VELOCITY = 3
INLET_ALPHA = 0.2

# exit pressure boundary condition:
EXIT_PRESSURE = P_ATM


inj_locs = np.array([NX//2])  # Inject at midpoint and 3/4 length
inj_rates = np.array([0.0001])   # kg/s of gas at each location


# %%
# optimization routine for implicit solves
# (quadratic) obj: min ||Ax-b||_2^2
fn_f = lambda A, x, b: x.T @ A.T @ A @ x - 2 * b.T @ A @ x + b @ b

# gradient and hessian
fn_g = lambda A, x, b: 2 * A.T @ (A @ x - b)
fn_h = lambda A, x, b: 2 * A.T @ A

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

        step = lr * np.linalg.inv(h) @ g  # newton step (for minimizing)
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

    
# %%

def fn_fric(rho, v, mu):
    """Calculates Blasius friction factor."""
    Re_m = rho * v * D / mu
    fB = 0.316 / (np.abs(Re_m+0.001)**0.25)
    # print(Re_m, rho, v, D, mu)
    return fB * np.sign(Re_m)

def calc_u_star_steady(p, rho, mu, u_guess, verbose:bool=False) -> dict:
    """
    Docstring for calc_u_star_steady
    
    :param p: Description
    :param rho: Description
    :param mu: Description
    :param verbose: Description
    :type verbose: bool
    :return: Description
    :rtype: dict[Any, Any]
    """
    def system(x):
        # self convection term: upwind
        # pressure term: center difference
        A = np.zeros((NX, NX))
        b = np.zeros((NX))
        for i in range(NX):
            if i == 0:
                # for inlet velocity boundary
                A[i, i] = 1
                b[i] = INLET_VELOCITY
            elif i == NX - 1:
                # for exit pressure boundary
                # enforce no velocity gradient
                # pressure value is handled in pressure corrector step
                A[i, i] = 1
                A[i, i - 1] = - 1
                b[i] = 0
            else:
                A[i, i] = x[i] / dx
                A[i, i - 1] = - x[i] / dx
                b[i] = - (1/rho[i]) * (p[i + 1] - p[i - 1]) / (2 * dx) - g_eff - fn_fric(rho[i], x[i], mu[i])/(2*D) * x[i]**2
        return A, b
    x = min_system(system, x0=u_guess, verbose=verbose)
    return x

def calc_u_star_transient(x, p, rho, mu) -> dict:
    """
    Docstring for calc_u_star_transient
    
    :param x: Description
    :param p: Description
    :param rho: Description
    :param mu: Description
    :return: Description
    :rtype: dict[Any, Any]
    """
    x_next = np.zeros(NX)
    for i in range(NX):
        if i == 0:
            x_next[i] = INLET_VELOCITY
        elif i == NX - 1:
            x_next[i] = x_next[i-1]
        else:
            RHS = - (1/rho[i]) * (p[i+1] - p[i-1])/(2*dx) - g_eff - fn_fric(rho[i], x[i], mu[i])/(2*D) * x[i]**2
            x_next[i] = x[i] + dt * (-x[i] * (x[i] - x[i - 1]) / dx - RHS)
    return {"x": x_next}

# %%
def calc_pressure_correction(u, rho, p, rlx=dt, rho_last=None, injection_locations=None, injection_rates=None) -> dict:
    """"""
    def system(x):
        A = np.zeros((NX, NX))
        b = np.zeros((NX))
        for i in range(NX):
            if i == 0:  # no pressure gradient at inlet - velocity boundary
                A[i, i] = 1
                A[i, i + 1] = -1
                b[i] = 0
            elif i == NX - 1:  # no correction at exit - pressure boundary
                A[i, i] = 1
                b[i] = 0
            else:  # center difference for interior
                A[i, i + 1] = (rlx / dx**2)
                A[i, i] = (rlx / dx**2) * (-2)
                A[i, i - 1] = (rlx / dx**2)

                if False:
                    b[i] = (rho[i + 1] * u[i + 1] - rho[i - 1] * u[i - 1]) / (2 * dx)
                    #b[i] = (rho[i] * u[i] - rho[i - 1] * u[i - 1]) / dx
                    if rho_last is not None:  # include transient term on RHS
                        b[i] += (rho[i] - rho_last[i]) / dt
                else:
                    # Rhie-Chow interpolation for face velocities
                    # Coefficient at faces (averaged)
                    coeff_plus = 0.5 * (dx / (rho[i] * np.abs(u[i])) + dx / (rho[i+1] * np.abs(u[i+1])))
                    coeff_minus = 0.5 * (dx / (rho[i-1] * np.abs(u[i-1])) + dx / (rho[i] * np.abs(u[i])))
                    
                    # Face velocity at i+1/2 with Rhie-Chow correction
                    u_interp_plus = 0.5 * (u[i] + u[i+1])
                    dp_face_plus = (p[i+1] - p[i]) / dx
                    dp_cell_avg_plus = 0.5 * ((p[i+1] - p[i-1])/(2*dx) + (p[i+2] - p[i])/(2*dx)) if i < NX-2 else (p[i+1] - p[i-1])/(2*dx)
                    u_face_plus = u_interp_plus - coeff_plus * (dp_face_plus - dp_cell_avg_plus)
                    
                    # Face velocity at i-1/2 with Rhie-Chow correction
                    u_interp_minus = 0.5 * (u[i-1] + u[i])
                    dp_face_minus = (p[i] - p[i-1]) / dx
                    dp_cell_avg_minus = 0.5 * ((p[i] - p[i-2])/(2*dx) + (p[i+1] - p[i-1])/(2*dx)) if i > 1 else (p[i+1] - p[i-1])/(2*dx)
                    u_face_minus = u_interp_minus - coeff_minus * (dp_face_minus - dp_cell_avg_minus)
                    
                    # Density at faces (averaged)
                    rho_face_plus = 0.5 * (rho[i] + rho[i+1])
                    rho_face_minus = 0.5 * (rho[i-1] + rho[i])
                    
                    # Mass flux divergence using face velocities
                    b[i] = (rho_face_plus * u_face_plus - rho_face_minus * u_face_minus) / dx
                    
                    if rho_last is not None:  # include transient term on RHS
                        b[i] += (rho[i] - rho_last[i]) / dt

                if injection_locations is not None and i in injection_locations and injection_rates is not None:
                    idx = np.where(injection_locations == i)[0][0]
                    m_dot_inj = injection_rates[idx]  # kg/s injected
                    A_pipe = np.pi * (D/2)**2
                    source = m_dot_inj / (A_pipe * dx)
                    b[i] += source

        return A, b
    A, b = system(None)
    p_corr = np.linalg.solve(A, b)
    return {"x": p_corr}
    return min_system(system, x0=alpha, verbose=False)

# %%
def void_transport(u_g, rho_g, alpha, alpha_last=None, injection_locations=None, injection_rates=None):
    def system(x):
        A = np.zeros((NX, NX))
        b = np.zeros((NX))
        for i in range(NX):
            if i == 0:  # fixed void at inlet
                A[i, i] = 1
                b[i] = INLET_ALPHA
            elif i == NX - 1:  # no void gradient at exit
                A[i, i] = 1
                A[i, i - 1] = - 1
                b[i] = 0
            else:
                # A[i, i + 1] = rho_g[i + 1] * u_g[i + 1] / (2 * dx)
                # A[i, i - 1] = - rho_g[i - 1] * u_g[i - 1] / (2 * dx)
                A[i, i] = rho_g[i] * u_g[i] / dx
                A[i, i - 1] = - rho_g[i - 1] * u_g[i - 1] / dx
                b[i] = 0
                # Source term for gas injection
                if injection_locations is not None and i in injection_locations and injection_rates is not None:
                    idx = np.where(injection_locations == i)[0][0]
                    m_dot_inj = injection_rates[idx]  # kg/s injected
                    A_pipe = np.pi * (D/2)**2
                    source = m_dot_inj / (A_pipe * dx)
                    b[i] += source
                if alpha_last is not None:
                    b[i] += (alpha_last[i] - alpha[i]) / dt
        return A, b
    A, b = system(None)
    #sol = np.linalg.inv(A.T @ A) @ A.T @ b
    sol = np.linalg.solve(A, b)
    return {"x": sol}
    return min_system(system, x0=alpha, verbose=False)

# %%
# variables
alpha = np.full(NX, INLET_ALPHA)
u_m = np.full(NX, INLET_VELOCITY)
u_g = u_g = C0 * u_m + V_gj
u_l = u_l = (u_m - alpha * u_g) / (1 - alpha)
p = np.full(NX, float(P_ATM))
rho_g = np.full(NX, fn_rho_g(P_ATM, T_AMB))

# convergence
RLX_P = 0.1
RLX_U = 0.2

# set P to single phase liquid hydrostatic head
p += (X[-1] - X) * g_eff * rho_l

#print(p, rho_g)

rho_m = alpha*rho_g + (1-alpha)*rho_l
mu_m = alpha*mu_g + (1-alpha)*mu_l

p_corr_log = []

# walkthrough of simple loop
for itr in range(1000):
    u_star_sol = calc_u_star_steady(p, rho_m, mu_m, u_guess=u_m, verbose=False)
    u_star = u_star_sol["x"]
    u_star[0] = INLET_VELOCITY
    if "converged" in u_star_sol.keys():
        if not u_star_sol["converged"]:
            print(f"Warning at iteration {itr}: u_star not converged.")
    if False:
        plt.figure(figsize=(6,2))
        plt.title(f"Residual: u*\nConverged: {u_star_sol['converged']}")
        plt.semilogy(u_star_sol["r_log"])

    # calculate pressure correction
    p_corr = calc_pressure_correction(u_star, rho_m, p)# injection_locations=inj_locs, injection_rates=inj_rates)
    if False:
        fig, ax1 = plt.subplots(figsize=(6,2))
        ax2 = ax1.twinx()
        plt.title("Pressure")
        ax1.plot(p)
        ax2.plot(p_corr["x"])
    p += p_corr["x"] * RLX_P
    if "converged" in p_corr.keys():
        if not p_corr["converged"]:
            print(f"Warning at iteration {itr}: p_corr not converged.")
    p[-1] = P_ATM
    p_corr_log.append(np.linalg.norm(p_corr["x"]))

    u_corr = np.zeros_like(u_star)
    # u_corr[1:-1] = - 1 / rho_m[1:-1] * (p_corr["x"][2:] - p_corr["x"][:-2]) / (2 * dx)  # transient correction, oops
    u_corr[1:-1] = - dx / (rho_m[1:-1] * INLET_VELOCITY) * (p_corr["x"][2:] - p_corr["x"][:-2]) / (2 * dx)  # random coefficients go!
    u_m = u_star + u_corr * RLX_U
    u_m[0] = INLET_VELOCITY

    # update gas velocity
    u_g = C0 * u_m + V_gj

    # solve void transport using known gas velocity
    void_update_solution = void_transport(u_g, rho_g, alpha, injection_locations=inj_locs, injection_rates=inj_rates)
    if "converged" in void_update_solution.keys():
        if not void_update_solution["converged"]:
            print(f"Warning at iteration {itr}: void_update_solution not converged.")
    alpha = void_update_solution["x"]
    alpha[0] = INLET_ALPHA

    alpha = np.clip(alpha, a_min=0.01, a_max=0.99)
    #plt.figure()
    #plt.title("alpha debug")
    #plt.plot(alpha)

    # update liquid velocity
    u_l = u_m - alpha * u_g
    u_l /= (1 - alpha)

    # update gas density via pressure field
    rho_g = fn_rho_g(p, T_AMB)

    # update mixture properties
    rho_m = alpha*rho_g + (1-alpha)*rho_l
    mu_m = alpha*mu_g + (1-alpha)*mu_l

plt.figure()
plt.semilogy(p_corr_log)
plt.title("Residuals")
plt.xlabel("Iteration")
plt.ylabel("Pressure Correction Magnitude")

# %%
plt.figure(figsize=(6,3))
plt.plot(X, u_m, label="Mixture Velocity")
plt.plot(X, u_g, label="Gas Velocity")
plt.plot(X, u_l, label="Liquid Velocity")
plt.legend()
plt.xlabel("Axial Location (m)")
plt.ylabel("Velocity (m/s)")

# %%
fig, ax1 = plt.subplots(figsize=(6,3))
ax2 = ax1.twinx()
ax1.plot(X, rho_g, label="gas density")
ax2.plot(X, alpha, color="red", linestyle="--", label="void fraction")
fig.legend()
ax1.set_xlabel("Axial Location (m)")
ax1.set_ylabel("Density (kg/m^3)")

# %%
plt.figure(figsize=(6,3))
plt.plot(X, p/1e3, label="Pressure")
plt.plot(X, ((X[-1] - X) * g_eff * rho_l + P_ATM)/1e3, label="Hydrostatic Head")
plt.axhline(P_ATM/1e3, color="gray", linestyle="--", label="Atmospheric Pressure")
plt.legend()
plt.xlabel("Axial Location (m)")
plt.ylabel("Pressure (kPa)")

# %%
fig, ax1 = plt.subplots(figsize=(6,3))
ax2 = ax1.twinx()
ax1.plot(rho_g*alpha*u_g + rho_l*(1-alpha)*u_l, color="purple")
ax2.plot(rho_g*alpha + rho_l*(1-alpha), color="blue")
ax2.axhline(rho_l, color="red")

# %%
A_pipe = np.pi * (D/2)**2

# Mass flow rates (kg/s)

inlet_mass_flow = (rho_g[0] * u_g[0] * alpha[0] + rho_l * u_l[0] * (1 - alpha[0])) * A_pipe
outlet_mass_flow = (rho_g[-1] * u_g[-1] * alpha[-1] + rho_l * u_l[-1] * (1 - alpha[-1])) * A_pipe
total_injected = np.sum(inj_rates)

expected_outlet = inlet_mass_flow + total_injected
error = np.abs(outlet_mass_flow - expected_outlet) / expected_outlet

print(f"Inlet mass flow: {inlet_mass_flow:.6f} kg/s")
print(f"Injected mass: {total_injected:.6f} kg/s")
print(f"Expected outlet: {expected_outlet:.6f} kg/s")
print(f"Actual outlet: {outlet_mass_flow:.6f} kg/s")
print(f"Conservation error: {error*100:.3f}%")
# %%
