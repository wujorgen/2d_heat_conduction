# Some 1D Heat Conduction Stuff
## 1D Steady-State Heat Conduction -> steady1d.c
assumes node 0 has a heat flux boundary condition
assumes node N-1 has a convection boundary condition
assumes inner nodes only conduct to adjacent nodes and do not convect, radiate
## 1D Transient Heat Conduction -> fullyimplicit1d.c
same assumptions as the steady heat conduction, just transient.
also uses a fully implicit discretization scheme.
### implicit vs explicit
explicit discretization scheme: a point is fully defined by the known surrounding temperature from the previous timestep.
implicit discretization scheme: a point is defined by only its previous temp and the other unknown temps in the current timestep.
The explicit discretization is obviously convienent but is subject to stability criterion among other issues. fully implicit >>>