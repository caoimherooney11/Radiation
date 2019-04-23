from firedrake import *
import numpy as np
import matplotlib.pyplot as plt
from radiation_cell_problem import solveCellProblem
from effective_conductivity import effectiveConductivity
from interpolate_keff import keff


def generate_keff(mesh, radius, k):
    length = 10 # set reasonable range of T values
    lists = [ [] for _ in range(5) ] 
    out = File("Output/Psi.pvd")
    for T in np.linspace(0, 1, length):
        eps = 1e-10
        T = 1. * (T + eps)
        Psi = solveCellProblem(mesh, radius, k, T)
        #Psi = out.write(Psi, time = T)

        (rad_integral, cond_integral) = effectiveConductivity(Psi, k, T)
        lists[0].append(T)
        l = 1
        for i in range(2):
            for j in range(2):
                lists[l].append(rad_integral[i][j] + cond_integral[i][j])
                l = l + 1
    
    lists_invert = list(map(list, zip(*lists)))
    np.savetxt("Output/BigDatasets/effective_conductivity.csv", lists_invert, delimiter=",")
    
    plt.plot(lists[0], lists[1], 'ro', lists[0], keff(lists[0], lists, 1), 'y')
    plt.savefig("Output/effective_conductivity.png", bbox_inches="tight")

    return lists



