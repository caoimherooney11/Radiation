from firedrake import *
import numpy as np
import matplotlib.pyplot as plt
from radiation_cell_problem import solveCellProblem
from effective_conductivity import effectiveConductivity
from interpolate_keff import keff


def generate_keff(data_name, mesh, radius, k, tau, c, vf, nonlinear, dim, dimensional=False):
    if nonlinear:
        length = 10 # set reasonable range of T values
        if dimensional:
            limit = 6500 # set highest T to consider (K)
            #limit = 4 # set highest T to consider (K)
        else:
            limit = 4 # set highest T to consider
        lists = [ [] for _ in range(5) ] 
        out = File("Output/Psi.pvd")
        for T in np.linspace(0, limit, length):
            warning("solving for T = %f" % T)
            eps = 1e-10
            T = 1. * (T + eps)
            Psi = solveCellProblem(mesh, radius, k, T, tau, c, vf, nonlinear, dim)
            #Psi = out.write(Psi, time = T)

            (rad_integral, cond_integral) = effectiveConductivity(Psi, k, T, radius, tau, c, vf, nonlinear)
            lists[0].append(T)
            l = 1
            for i in range(2):
                for j in range(2):
                    val = rad_integral[i][j] + cond_integral[i][j]
                    if dimensional: 
                        val = val * k
                    lists[l].append(val)
                    l = l + 1
        
        lists_invert = list(map(list, zip(*lists)))
        np.savetxt(data_name + ".csv", lists_invert, delimiter=",")
        
        plt.plot(lists[0], lists[1], 'ro', lists[0], keff(lists[0], lists, limit, 1), 'y')
        plt.savefig(data_name + ".png", bbox_inches="tight")

        return lists

    else:
        Psi = solveCellProblem(mesh, radius, k, None, None, c, vf, nonlinear, dim)
        (rad_integral, cond_integral) = effectiveConductivity(Psi, k, None, radius, None, c, vf,  nonlinear)
        integral = np.zeros((2,2))
        for i in range(2):
            for j in range(2):
               integral[i, j] = rad_integral[i][j] + cond_integral[i][j]
               if dimensional:
                   integral[i, j] = integral[i, j] * k
        
        return integral




