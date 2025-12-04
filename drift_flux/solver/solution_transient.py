import matplotlib.pyplot as plt
import numpy as np

from classes import Solution, CollocatedSolution
from constants import P_ATM, T_AMB, T_end, dt, NT, dx, N_CELLS, N_FACES, CELLS, FACES, g_eff, rho_l, fn_rho_g, mu_g, mu_l, C0, fn_Vgj, L
from functions import calc_U_star_transient, calc_pressure_correction_transient, gas_continuity_for_void_transient
from helper_functions import fn_check_step_convergence


def transient_solution(initial_condition:Solution|None=None):
    """Solves time-varying mixture continuity and momentum equations using the SIMPLE algorithm.
        A pressure perturbation at the inlet is applied from 0.1 to 0.5 seconds.
    """

    # inlet velocity boundary condition
    INLET_VELOCITY = 2
    INLET_ALPHA = 0.2

    # exit pressure boundary condition:
    INLET_PRESSURE = P_ATM + rho_l * 9.81 * L
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

    # previous timestep variables for transient loop
    rho_last = rho_m
    rho_g_last = rho_g
    alpha_last = alpha
    u_last = u_m

    # convergence
    RLX_P = 0.1
    RLX_U = 0.4
    MAX_ITR = 500

    CONVERGED = False
    CTOL = 5e-3

    corr_log = []
    solution_log = []

    for tdx in range(1):
        # reset convergence tracker for time step
        CONVERGED = False
        # prepare convergence log
        corr_log.append([])

        for itr in range(MAX_ITR):
            u_star_sol = calc_U_star_transient(p, rho_m, mu_m, u_m)
            u_star = u_star_sol["x"]
            # enforce no velocity gradient at inlet and exit
            u_star[0] = u_star[1]
            u_star[-1] = u_star[-2]

            # calc pressure correction
            p_corr_sol = calc_pressure_correction_transient(u=u_star, rho=rho_m, rho_last=rho_last)
            fn_check_step_convergence(p_corr_sol, itr, "p_corr_sol")
            p_corr = p_corr_sol["x"]
            corr_log[tdx].append(np.linalg.norm(p_corr))

            if np.linalg.norm(p_corr) < CTOL:
                CONVERGED = True
                break

            # calc velocity correction
            u_corr = -2 * dt / (rho_m[1:] + rho_m[:-1]) * (p_corr[1:] - p_corr[:-1])/dx

            # apply (and relax) corrections
            u_m = u_star + u_corr * RLX_U
            p += p_corr * RLX_P

            # enforce bc
            #u_m[0] = INLET_VELOCITY
            #p[-1] = 2 * P_ATM - p[-2]
            p[0] = INLET_PRESSURE
            p[1] = INLET_PRESSURE
            p[-1] = EXIT_PRESSURE
            p[-2] = EXIT_PRESSURE

            # update gas velocity via drift flux closure
            V_gj = fn_Vgj( (alpha[1:] + alpha[:-1])/2 )
            u_g = C0 * u_m + V_gj

            # solve gas continuity to update void fraction
            alpha_sol = gas_continuity_for_void_transient(rho_g=rho_g, u_g=u_g, void_in=INLET_ALPHA, void_last=alpha_last, rho_g_last=rho_g_last)
            fn_check_step_convergence(alpha_sol, itr, "alpha_sol")
            alpha = alpha_sol["x"]
            #alpha[0] = INLET_ALPHA  # wrong - handled in gas continuity solve
            alpha = np.clip(alpha, a_min=0.01, a_max=0.99)  # TODO: diff values?

            # update liquid velocity
            u_l = u_m - (alpha[1:] + alpha[:-1])/2 * u_g
            u_l /= (1 - (alpha[1:] + alpha[:-1])/2)

            # update gas density via pressure field
            rho_g = fn_rho_g(p, T_AMB)

            # update mixture properties
            rho_m = alpha*rho_g + (1-alpha)*rho_l
            mu_m = alpha*mu_g + (1-alpha)*mu_l

        print("after timestep loop - converged:", CONVERGED)
        plt.figure()
        plt.semilogy(corr_log[tdx])
        plt.show()
        breakpoint()
        # after timestep converges, update 
        rho_last = rho_m
        rho_g_last = rho_g
        alpha_last = alpha
        u_last = u_m


if __name__ == "__main__":
    print("NOT FINISHED YET")
    # TODO: this should generate a steady state solution to feed into the transient. probably.
    transient_solution()
