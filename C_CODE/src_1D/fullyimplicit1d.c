#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Functions.h"
#include "Node.h"

int main()
{
	int N = 11;
	struct Node *Mesh;
	Mesh = (struct Node *)calloc(N, sizeof(struct Node));
	// BC_Mesh
	if (Mesh == NULL)
	{
		printf("FUCK NAH\n");
		return -1;
	}
	double T_0[N]; // fill initial condition
	for (int i = 0; i < N; i++)
	{
		T_0[i] = 20 + (double)(i);
	}
	// Fill thermal conductivities and step sizes
	// k_steel = 45 Watts / (Kelvin Meter)
	// delta_x = 0.1 Meters i guess
	// rho_steel = 7850 kg/m^3
	double rho = 7850;
	// c = 420 J / (kg Kelvin) -> specific heat capacity
	double specific_heat = 420;
	// Set initial guess for all according to temperature initial condition
	for (int i = 0; i < N; i++)
	{
		Mesh[i].T = T_0[i];
		Mesh[i].k = 45.0;
		Mesh[i].dx = 0.05;
	}

	// Heat flux boundary condition
	// strcpy(Mesh[0].bc, "HeatFlux");
	double phi_bound = 1000; // Watts / Square Meter
	// Convection boundary condition
	// strcpy(Mesh[N - 1].bc, "Convection");
	double T_amb = 22;
	double h_amb = 25.32; // convection, steel-to-air, W / (m^2 K)

	// Print the initial guess temps
	printf("INITIAL: ");
	for (int i = 0; i < N; i++)
	{
		printf("%9.3f", Mesh[i].T);
	}
	printf("\n");

	// so the loop heirarchy here will be:
	// timestep -> iteration -> point-by-point

	// transient loop control variables
	double dt = 1.0; // timestep, seconds
	double transient_end_time = 1000000.0;
	// T_0 will be updated to show the temperatures at previous time step

	// convergence loop control variables
	double T_last[N]; // T_last shows previous iteration
	double stol = 0.001;
	bool ITERATE = true;
	double R;

	// equation coefficient declarations
	double a_p;
	double a_p0;
	double a_w;
	double a_e;

	// TRANSIENT CONTROL LOOP
	for (double t = 0; t <= transient_end_time; t += dt)
	{
		// store temps from last timestep.
		// redundant for iteration zero but eh
		for (int i = 0; i < N; i++)
		{
			T_0[i] = Mesh[i].T;
		}
		// ITERATION CONTROL LOOP
		ITERATE = true;
		for (int iter = 0; ITERATE; iter++)
		{
			// store temps from last iteration
			for (int i = 0; i < N; i++)
			{
				T_last[i] = Mesh[i].T;
			}
			// visit mesh nodes: POINT-BY-POINT
			for (int i = 0; i < N; i++)
			{
				a_p0 = rho * specific_heat * Mesh[i].dx / dt;
				if (i == 0)
				{
					a_w = 0.0;
					a_e = h_mean(Mesh[i].k, Mesh[i + 1].k) / Mesh[i].dx;
					a_p = a_w + a_e + a_p0;
					Mesh[i].T = a_e * Mesh[i + 1].T + phi_bound + a_p0 * T_0[i];
					Mesh[i].T /= a_p;
				}
				else if (i == N - 1)
				{

					a_w = h_mean(Mesh[i - 1].k, Mesh[i].k) / Mesh[i].dx;
					a_e = 0.0;
					a_p = a_w + a_e + a_p0 + h_amb;
					Mesh[i].T = a_w * Mesh[i - 1].T + h_amb * T_amb + a_p0 * T_0[i];
					Mesh[i].T /= a_p;
				}
				else
				{
					a_w = h_mean(Mesh[i - 1].k, Mesh[i].k) / Mesh[i].dx;
					a_e = h_mean(Mesh[i].k, Mesh[i + 1].k) / Mesh[i].dx;
					a_p = a_w + a_e + a_p0;
					Mesh[i].T = a_w * Mesh[i - 1].T + a_e * Mesh[i + 1].T + a_p0 * T_0[i];
					Mesh[i].T /= a_p;
				}
			}
			// check for convergence - I think this is basically L_inf
			ITERATE = false;
			for (int i = 0; i < N; i++)
			{
				if (fabs(Mesh[i].T - T_last[i]) > stol)
				{
					ITERATE = true;
				}
			}
			// print temps after point-by-point, every nth iteration OR when converged
			if (!ITERATE) // iter % 1 == 0 ||
			{
				//if (!ITERATE)
				//	printf("CONVERGED - SOLUTION BELOW FOR TIMESTEP %f:\n", t + dt);
				printf("time %5.2f ITR %d: ", t + dt, iter);
				for (int i = 0; i < N; i++)
				{
					printf("%9.3f", Mesh[i].T);
				}
				printf("\n");
			}
		}
	}
	return 0;
}
