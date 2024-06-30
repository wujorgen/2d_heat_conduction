#ifndef NODE_H
#define NODE_H

struct MegaNode
{
    double T; // Temperature
    double k; // Thermal Conductivity
    double dx; // Node Width
    char bc[10]; // "None", "Temp", "HeatFlux", "Convection"
    double T_bc; // Temperature Boundary Condition
    double q_bc; // Heat Flux Bondary Condition
    double h_bc; // Convection Heat Transfer Coefficient Boundary Condition
};

struct Node
{
    double T;
    double k;
    double a;
    double delta_x;
};


#endif