#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Functions.h"
#include "Node.h"

int main()
{
    int N = 9;
    struct MegaNode *Mesh;
    Mesh = (struct MegaNode *)calloc(N, sizeof(struct MegaNode));
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
    for (int i = 0; i < N; i++)
    {
        Mesh[i].T = 20;
        Mesh[i].k = 45.0;
        Mesh[i].dx = 0.1;
        strcpy(Mesh[i].bc, "None");
    }

    // Heat flux boundary condition
    strcpy(Mesh[0].bc, "HeatFlux");
    double phi_bound = 100; // Watts / Square Meter
    // Convection boundary condition
    strcpy(Mesh[N - 1].bc, "Convection");
    double T_amb = 22;
    double h_amb = 25.32; // convection, steel-to-air, W / (m^2 K)

    // Print the initial guess temps
    printf("INITIAL: ");
    for (int i = 0; i < N; i++)
    {
        printf("%9.3f", Mesh[i].T);
    }
    printf("\n");

    return 0;
}