#include <stdio.h>
#include <stdlib.h>

#include "Functions.h"

double h_mean(double a, double b) // harmonic mean
{
    return 2.0 * a * b / (a + b);
    //return (a + b) / 2.0;
}