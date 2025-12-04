import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass


@dataclass
class Solution:
    N_FACES: int
    N_CELLS: int
    MESH: NDArray[np.float64]
    U: NDArray[np.float64]
    U_L: NDArray[np.float64]
    U_G: NDArray[np.float64]
    RHO: NDArray[np.float64]
    RHO_G: NDArray[np.float64]
    P: NDArray[np.float64]
    CONVERGED: bool
    LOG: NDArray[np.float64]

    def cells_to_faces(self, field:str) -> NDArray[np.float64]|None:
        try:
            arr = self.__getattribute__(str(field).upper())
            arr = (arr[1:] + arr[:-1]) / 2
        except AttributeError as e:
            print(e)
            arr = None
        finally:
            return arr


class CollocatedSolution:
    def __init__(self, sol:Solution) -> None:
        self.NX: int = sol.N_FACES
        self.X: NDArray[np.float64] = sol.MESH
        self.U: NDArray[np.float64] = sol.U
        self.U_L: NDArray[np.float64] = sol.U_L
        self.U_G: NDArray[np.float64] = sol.U_G
        self.RHO: NDArray[np.float64]|None = sol.cells_to_faces("RHO")
        self.RHO_G: NDArray[np.float64]|None = sol.cells_to_faces("RHO_G")
        self.P: NDArray[np.float64]|None = sol.cells_to_faces("P")

    def plot(self, variable:str, show:bool=False, save:str|None=None) -> None:
        try:
            arr = self.__getattribute__(str(variable).upper())
        except AttributeError as e:
            print(e)
            return
        plt.figure()
        plt.plot(self.X, arr)
        if show:
            plt.show()
        if isinstance(save, str):
            plt.savefig(save)
        plt.close()
