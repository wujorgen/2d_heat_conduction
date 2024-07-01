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
    double h_amb; // Convection Heat Transfer Coefficient Boundary Condition
    double T_amb; // Ambient Temperature
};

struct Node
{
    double T;
    double k;
    double a;
    double dx;
};

struct SuperNode
{
    double T; // Temperature
    double k; // Thermal Conductivity
    double dx; // Node Width, X
    double dy; // Node Width, Y
    // Boundary Condition Types: "None", "Temp", "HeatFlux", "Convection"
    char bc_W[10]; // West / Left
    char bc_E[10]; // East / Right
    char bc_N[10]; // North / Up
    char bc_S[10]; // South / Down
    // Boundary Condition Constants:

};

#endif