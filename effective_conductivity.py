from firedrake import *
import numpy as np

def effectiveConductivity(Psi, k, T, radius, tau, c, VF, nonlinear):
    mesh = Psi.ufl_domain()
    length = assemble(Constant(1.0) * ds(1, domain=mesh))
    #vf = 1 / length
    vf = VF(radius)
    area = assemble(Constant(1.0) * dx(domain=mesh))
    dim = len(Psi)
    X = SpatialCoordinate(mesh)

    def G(u):
        return u - vf * assemble(u * ds(1, domain=mesh))

    rad_integral = np.zeros((dim, dim))
    cond_integral = np.zeros((dim, dim))

    if nonlinear:
        alpha = 4 * c * (T + tau)**3
    else:
        alpha = c
    for i in range(dim):
        for j in range(dim):
            rad_integral[i,j] = alpha * assemble(X[i] * G(Psi[j] + X[j]) * ds(1, domain=mesh)) / area
            #cond_integral[i,j] = k * assemble((Identity(2)[i,j] + grad(Psi)[i,j]) * dx) / area
            cond_integral[i,j] = assemble((Identity(2)[i,j] + grad(Psi)[i,j]) * dx) / area

    return (rad_integral, cond_integral)    
