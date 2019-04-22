from firedrake import *
import numpy as np

def effectiveConductivity(Psi, k, T):
    sigma = 1.
    mesh = Psi.ufl_domain()
    area = assemble(Constant(1.0) * ds(1, domain=mesh))
    vf = 1 / area
    dim = len(Psi)
    X = SpatialCoordinate(mesh)

    def G(u):
        return u - vf * assemble(u * ds(1, domain=mesh))

    rad_integral = np.zeros((dim, dim))
    cond_integral = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            rad_integral[i,j] = 4 * sigma * T**3 * assemble(X[i] * G(Psi[j] + X[j]) * ds(1, domain=mesh))
            cond_integral[i,j] = k * assemble((Identity(2)[i,j] + grad(Psi)[i,j]) * dx)

    return (rad_integral, cond_integral)    
