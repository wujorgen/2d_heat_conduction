#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Functions.h"
// #include "Mesh.h"
#include "Node.h"

int main()
{
	int N = 9;
	struct Node *Mesh;
	Mesh = (struct Node *)calloc(N, sizeof(struct Node));
	// BC_Mesh
	if (Mesh == NULL)
	{
		printf("FUCK NAH\n");
		return -1;
	}
	// Fill thermal conductivities and step sizes
	// k_steel = 45 Watts / (Kelvin Meter)
	// delta_x = 0.1 Meters i guess
	// Set boundary condition types to None
	// Set initial guess for all temperatures to 0 (Celsius I think)
	// Temperature initial guess doesn't really matter here since its steady state
	for (int i = 0; i < N; i++)
	{
		Mesh[i].T = 20;
		Mesh[i].k = 45.0;
		Mesh[i].dx = 0.1;
		// strcpy(Mesh[i].bc, "None");
	}

	// Heat flux boundary condition
	// strcpy(Mesh[0].bc, "HeatFlux");
	double phi_bound = 100; // Watts / Square Meter
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

	// Point-by-point
	// Temperatures from previous iteration
	double T_last[N];

	// solution tolerance
	double stol = 0.001;
	bool ITERATE = true;

	// residual tracker, L_inf or L2
	double R;

	// equation coefficient declarations
	double a_p;
	double a_w;
	double a_e;

	for (int iter = 0; ITERATE; iter++) // iterate
	{
		// store temps from last iteration
		for (int i = 0; i < N; i++)
		{
			T_last[i] = Mesh[i].T;
		}
		// visit mesh point by point
		for (int i = 0; i < N; i++)
		{
			if (i == 0)
			{
				a_w = 0.0;
				a_e = h_mean(Mesh[i].k, Mesh[i + 1].k) / Mesh[i].dx;
				a_p = a_w + a_e;
				Mesh[0].T = a_e * Mesh[1].T + phi_bound;
				Mesh[0].T /= a_p;
			}
			else if (i == N - 1)
			{
				a_w = h_mean(Mesh[i - 1].k, Mesh[i].k) / Mesh[i].dx;
				a_e = 0.0;
				a_p = a_w + a_e + h_amb;
				Mesh[i].T = a_w * Mesh[i - 1].T + h_amb * T_amb;
				Mesh[i].T /= a_p; //(h_mean(Mesh[0].k, Mesh[1].k) / Mesh[0].dx) + h_amb;
			}
			else
			{
				a_w = h_mean(Mesh[i - 1].k, Mesh[i].k) / Mesh[i].dx;
				a_e = h_mean(Mesh[i].k, Mesh[i + 1].k) / Mesh[i].dx;
				a_p = a_w + a_e;
				Mesh[i].T = a_w * Mesh[i - 1].T + a_e * Mesh[i + 1].T;
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
		// print temps after point-by-point, every tenth iteration OR when converged
		if (iter % 20 == 0 || !ITERATE)
		{
			if (!ITERATE)
				printf("CONVERGED - SOLUTION BELOW:\n");
			printf("ITR %d: ", iter);
			for (int i = 0; i < N; i++)
			{
				printf("%9.3f", Mesh[i].T);
			}
			printf("\n");
		}
	}

	return 0;
}
