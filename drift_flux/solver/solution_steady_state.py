import matplotlib.pyplot as plt
import numpy as np

from classes import Solution, CollocatedSolution
from constants import P_ATM, T_AMB, dx, N_CELLS, N_FACES, CELLS, FACES, g_eff, rho_l, fn_rho_g, mu_g, mu_l, C0, fn_Vgj
from functions import calc_U_star_steady, calc_pressure_correction_steady, gas_continuity_for_void_steady
from helper_functions import fn_check_step_convergence

def steady_state_solution():
    """Solves steady state mixture continuity and momentum equations using the SIMPLE algorithm."""
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

    for itr in range(MAX_ITR):
        # calc tentative velocity
        u_star_sol = calc_U_star_steady(p=p, rho=rho_m, mu=mu_m, u_last=u_m, u_in=INLET_VELOCITY)
        fn_check_step_convergence(u_star_sol, itr, "u_star_sol")
        u_star = u_star_sol["x"]
        u_star[0] = INLET_VELOCITY

        # calc pressure correction
        p_corr_sol = calc_pressure_correction_steady(u=u_star, rho=rho_m)
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
        alpha_sol = gas_continuity_for_void_steady(rho_g=rho_g, u_g=u_g, void_in=INLET_ALPHA, void_last=alpha)
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


if __name__ == "__main__":
    ss_sol = steady_state_solution()
    print(f"{ss_sol['converged']=}, {len(ss_sol['convergence_history'])}")

    plt.figure(figsize=(6,3))
    plt.semilogy(ss_sol["convergence_history"])
    plt.show()

    plt.figure(figsize=(6,3))
    plt.plot(ss_sol["u_m"], label="mixture")
    plt.plot(ss_sol["u_l"], label="liquid")
    plt.plot(ss_sol["u_g"], label="gas")
    plt.legend()
    plt.show()

    fig, ax1 = plt.subplots(figsize=(6,3))
    ax2 = ax1.twinx()
    ax1.plot(ss_sol["rho_g"], label="gas density")
    ax2.plot(ss_sol["alpha"], color="red", linestyle="--", label="void fraction")
    ax2.set_ylim(0,1)
    fig.legend()
    fig.tight_layout()
    plt.show()

    plt.figure(figsize=(6,3))
    plt.plot(ss_sol["p"], label="pressure")
    plt.plot((FACES[-1] - FACES) * g_eff * rho_l + P_ATM, label="hydrostatic pressure")
    plt.axhline(P_ATM, linestyle="--", color="green", label="Atmospheric")
    plt.legend()
    plt.show()


    inlet_mass_flux = ss_sol["rho_g"][0] * ss_sol["u_g"][0] * ss_sol["alpha"][0] + rho_l * ss_sol["u_l"][0] * (1 - ss_sol["alpha"][0])
    outlet_mass_flux = ss_sol["rho_g"][-1] * ss_sol["u_g"][-1] * ss_sol["alpha"][-1] + rho_l * ss_sol["u_l"][-1] * (1 - ss_sol["alpha"][-1])

    pct_mass_change = (outlet_mass_flux - inlet_mass_flux) / inlet_mass_flux * 100
    print(f"Difference in outlet to inlet mass flow: {pct_mass_change:.2f}%")

    mflux_in = ss_sol["rho_m"][0] * ss_sol["u_m"][0]
    mflux_out = ss_sol["rho_m"][-1] * ss_sol["u_m"][-1]
    print(f"Inlet: {mflux_in:.6f} kg/s/m^2, Outlet: {mflux_out:.6f} kg/s/m^2")
    print(f"Difference: {mflux_out - mflux_in:.6f} kg/s/m^2")
