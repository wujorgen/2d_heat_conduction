import numpy as np
from CoolProp.CoolProp import PropsSI


# Finite Volume Grid
N_FACES = 11
N_CELLS = N_FACES + 1  # include ghost cells
L = 1
dx = L / (N_FACES - 1)
FACES = np.linspace(0, L, N_FACES, endpoint=True)
CELLS = np.linspace(0 - dx / 2, L + dx / 2, N_CELLS, endpoint=True)

D = 1 / 39.37
A_pipe = np.pi * (D/2)**2
dt = 0.01
T_end = 1
NT = int(T_end / dt + 1)

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

# 1. Get the universal gas constant (J/mol/K)
R = PropsSI("gas_constant", "air")

# 2. Get the molar mass (kg/mol)
M = PropsSI("molemass", "air") # in kg/mol (Note: CoolProp uses kg/mol by default for this)

# 3. Calculate the specific gas constant (J/kg/K)
Rbar = R / M

fn_rho_g = lambda P, T: P / (Rbar * T)
mu_g = PropsSI("V", "P", P_ATM, "T", T_AMB, "Air")
rho_l = PropsSI("D", "P", P_ATM, "T", T_AMB, "Water")
mu_l = PropsSI("V", "P", P_ATM, "T", T_AMB, "Water")     


if __name__ == "__main__":
    print("Air Properties")
    print(f"Universal gas constant (R): {R:.3f} J/(mol·K)")
    print(f"Molar mass (M): {M:.3f} kg/mol")
    print(f"Specific gas constant (Rs): {Rbar:.3f} J/(kg·K)")
    print(f"Air density at room temp should be around 1.204 kg/m^3. Confirm: {fn_rho_g(P_ATM, T_AMB):.3f} kg/m^3")
    print(f"Water density at room temp: {rho_l:.3f} kg/m^3")
    print(f"Water dynamic viscosity at room temp: {mu_l:.3f} kg/m^3")
