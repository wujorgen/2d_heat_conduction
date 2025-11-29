# %%
import numpy as np
import scipy
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from dataclasses import dataclass
from copy import deepcopy

# Finite Volume Grid
N_FACES = 51
N_CELLS = N_FACES + 1  # include ghost cells
L = 1
dx = L / (N_FACES - 1)
FACES = np.linspace(0, L, N_FACES, endpoint=True)
CELLS = np.linspace(0 - dx / 2, L + dx / 2, N_CELLS, endpoint=True)

D = 1 / 39.37
A_pipe = np.pi * (D/2)**2
dt = 0.01
T_end = 1
NT = T_end / dt + 1

# orientation
theta = 90  # degrees
g_constant = 9.81
g_eff = g_constant * np.sin(np.deg2rad(theta))

# Drift flux parameters
C0 = 1.2
# V_gj = 0.35 * np.sqrt(g_eff * D)

fn_Vgj = lambda alpha: 0.35 * np.sqrt(g_eff * D) * (1 - alpha)**1  # exponent in [0.5, 1.5]

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

'''
# inlet velocity boundary condition
INLET_VELOCITY = 3
INLET_ALPHA = 0.2

# exit pressure boundary condition:
EXIT_PRESSURE = P_ATM

#inj_locs = np.array([NX//2])  # Inject at midpoint and 3/4 length
#inj_rates = np.array([0.0001])   # kg/s of gas at each location

# variables: cell centers
alpha = np.full(N_CELLS, INLET_ALPHA)

p = np.full(N_CELLS, float(P_ATM))
p += (CELLS[-1] - CELLS) * g_eff * rho_l  # set P to single phase liquid hydrostatic head

rho_g = fn_rho_g(p, T_AMB)
# rho_l = CONSTANT

rho_m = alpha*rho_g + (1-alpha)*rho_l
mu_m = alpha*mu_g + (1-alpha)*mu_l

# variables: cell faces
u_m = np.full(N_FACES, INLET_VELOCITY)
u_g = u_g = C0 * u_m + V_gj
u_l = u_l = (u_m - (alpha[1:] + alpha[:-1])/2 * u_g) / (1 - (alpha[1:] + alpha[:-1])/2)

# convergence
RLX_P = 0.1
RLX_U = 0.2
'''

# %%
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
                friction = fn_fric(rho_hat, x[i], mu_hat) / (2 * D) * x[i]**2
                b[i] = (- 1 / rho_hat) * (p[i+1] - p[i]) / dx - g_eff - friction
        return A, b
    sol = min_system(system, x0=u_last)
    return sol

def calc_U_star_transient(p, rho, mu, u_last):
    pass

# u_star_sol = calc_U_star_steady(p = p, rho = rho_m, mu = mu_m, u_last = u_m)

# %%
def calc_pressure_correction(u, rho, dt=dt, rho_last=None):
    if rho_last is None:
        dt = 1
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

                if rho_last is not None:
                    b[i] += (rho[i] - rho_last[i]) / dt
        return A, b
    #print(*system(None))
    sol = min_system(system, x0=np.zeros(N_CELLS))
    return sol
    return {"x": np.linalg.solve(*system(None))}

# p_corr_sol = calc_pressure_correction(u=u_star_sol["x"], rho=rho_m)

# %%
def gaa_continuity_for_void(rho_g, u_g, void_in, void_last, rho_g_last=None, injection_locations=None, injection_rates=None):
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
                # Source term for gas injection
                if (injection_locations is not None) and (i in injection_locations) and (injection_rates is not None):
                    idx = np.where(injection_locations == i)[0][0]
                    m_dot_inj = injection_rates[idx]  # kg/s injected
                    source = m_dot_inj / (A_pipe * dx)
                    b[i] += source
                if void_last is not None and rho_g_last is not None:
                    pass  # TODO: (rho_g*alpha - rho_g_last*alpha_last) / dt <- alpha needs to be moved to LHS, or assume alpha from last iteration
        return A, b
    sol = min_system(system, x0=void_last)
    return sol

# %%
def steady_state_solution():
    # inlet velocity boundary condition
    INLET_VELOCITY = 2
    INLET_ALPHA = 0.2

    # exit pressure boundary condition:
    EXIT_PRESSURE = P_ATM

    # variables: cell centers
    alpha = np.full(N_CELLS, INLET_ALPHA)

    p = np.full(N_CELLS, float(P_ATM))
    p += (CELLS[-1] - CELLS) * g_eff * rho_l  # set P to single phase liquid hydrostatic head

    rho_g = fn_rho_g(p, T_AMB)
    # rho_l = CONSTANT

    rho_m = alpha*rho_g + (1-alpha)*rho_l
    mu_m = alpha*mu_g + (1-alpha)*mu_l

    # variables: cell faces
    u_m = np.full(N_FACES, INLET_VELOCITY)
    u_g = u_g = C0 * u_m + fn_Vgj( (alpha[1:] + alpha[:-1])/2 )
    u_l = u_l = (u_m - (alpha[1:] + alpha[:-1])/2 * u_g) / (1 - (alpha[1:] + alpha[:-1])/2)

    # convergence
    RLX_P = 0.5
    RLX_U = 0.8
    MAX_ITR = 2000

    CONVERGED = False
    CTOL = 5e-3

    corr_log = []


    def fn_check_step_convergence(sol:dict, itr, name):
        if "converged" in sol.keys():
            if not sol["converged"]:
                print(f"Warning at iteration {itr}: {name} not converged.")


    for itr in range(MAX_ITR):
        # calc tentative velocity
        u_star_sol = calc_U_star_steady(p=p, rho=rho_m, mu=mu_m, u_last=u_m, u_in=INLET_VELOCITY)
        fn_check_step_convergence(u_star_sol, itr, "u_star_sol")
        u_star = u_star_sol["x"]
        u_star[0] = INLET_VELOCITY

        # calc pressure correction
        p_corr_sol = calc_pressure_correction(u=u_star, rho=rho_m)
        fn_check_step_convergence(p_corr_sol, itr, "p_corr_sol")
        p_corr = p_corr_sol["x"]
        corr_log.append(np.linalg.norm(p_corr))

        if np.linalg.norm(p_corr) < CTOL:
            CONVERGED = True
            break

        # calc velocity correction
        u_corr = -2 * dx / (rho_m[1:] + rho_m[:-1]) * (p_corr[1:] - p_corr[:-1])/dx

        # apply (and relax) corrections
        u_m = u_star + u_corr * RLX_U
        p += p_corr * RLX_P

        # enforce bc
        u_m[0] = INLET_VELOCITY
        #p[-1] = 2 * P_ATM - p[-2]
        p[-1] = P_ATM
        p[-2] = P_ATM

        # update gas velocity via drift flux closure
        V_gj = fn_Vgj( (alpha[1:] + alpha[:-1])/2 )
        u_g = C0 * u_m + V_gj

        # solve gas continuity to update void fraction
        alpha_sol = gaa_continuity_for_void(rho_g=rho_g, u_g=u_g, void_in=INLET_ALPHA, void_last=alpha)
        fn_check_step_convergence(alpha_sol, itr, "alpha_sol")
        alpha = alpha_sol["x"]
        # alpha[0] = INLET_ALPHA  # wrong - handled in gas continuity solve
        alpha = np.clip(alpha, a_min=0.01, a_max=0.99)  # TODO: diff values?

        # update liquid velocity
        u_l = u_m - (alpha[1:] + alpha[:-1])/2 * u_g
        u_l /= (1 - (alpha[1:] + alpha[:-1])/2)

        # update gas density via pressure field
        rho_g = fn_rho_g(p, T_AMB)

        # update mixture properties
        rho_m = alpha*rho_g + (1-alpha)*rho_l
        mu_m = alpha*mu_g + (1-alpha)*mu_l
    
    p_at_faces = (p[1:] + p[:-1]) / 2
    rho_g_at_faces = (rho_g[1:] + rho_g[:-1]) / 2
    rho_m_at_faces = (rho_m[1:] + rho_m[:-1]) / 2
    alpha_at_faces = (alpha[1:] + alpha[:-1]) / 2

    return {
        "u_m": u_m,
        "u_l": u_l,
        "u_g": u_g,
        "alpha": alpha_at_faces,
        "p": p_at_faces,
        "rho_m": rho_m_at_faces,
        "rho_g": rho_g_at_faces,
        "converged": CONVERGED,
        "convergence_history": corr_log,
    }

ss_sol = steady_state_solution()
print(f"{ss_sol['converged']=}, {len(ss_sol['convergence_history'])}")

# %%
plt.figure(figsize=(6,3))
plt.semilogy(ss_sol["convergence_history"])

plt.figure(figsize=(6,3))
plt.plot(ss_sol["u_m"], label="mixture")
plt.plot(ss_sol["u_l"], label="liquid")
plt.plot(ss_sol["u_g"], label="gas")
plt.legend()

fig, ax1 = plt.subplots(figsize=(6,3))
ax2 = ax1.twinx()
ax1.plot(ss_sol["rho_g"], label="gas density")
ax2.plot(ss_sol["alpha"], color="red", linestyle="--", label="void fraction")
ax2.set_ylim(0,1)
fig.legend()
fig.tight_layout()

plt.figure(figsize=(6,3))
plt.plot(ss_sol["p"], label="pressure")
plt.plot((FACES[-1] - FACES) * g_eff * rho_l + P_ATM, label="hydrostatic pressure")
plt.axhline(P_ATM, linestyle="--", color="green", label="Atmospheric")
plt.legend()

# %%
inlet_mass_flux = ss_sol["rho_g"][0] * ss_sol["u_g"][0] * ss_sol["alpha"][0] + rho_l * ss_sol["u_l"][0] * (1 - ss_sol["alpha"][0])
outlet_mass_flux = ss_sol["rho_g"][-1] * ss_sol["u_g"][-1] * ss_sol["alpha"][-1] + rho_l * ss_sol["u_l"][-1] * (1 - ss_sol["alpha"][-1])

pct_mass_change = (outlet_mass_flux - inlet_mass_flux) / inlet_mass_flux * 100
print(f"Difference in outlet to inlet mass flow: {pct_mass_change:.2f}%")

# %%
mflux_in = ss_sol["rho_m"][0] * ss_sol["u_m"][0]
mflux_out = ss_sol["rho_m"][-1] * ss_sol["u_m"][-1]
print(f"Inlet: {mflux_in:.6f} kg/s/m^2, Outlet: {mflux_out:.6f} kg/s/m^2")
print(f"Difference: {mflux_out - mflux_in:.6f} kg/s/m^2")

# %%
