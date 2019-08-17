from firedrake import *
import numpy as np
from interpolate_keff import keff

def solve_full(domain_dimensions, k_eff, BC, f):
    nx = domain_dimensions[0] 
    ny = domain_dimensions[1]
    n = 10
    mesh = RectangleMesh(200 * nx, 200 * ny, nx, ny)
    x = SpatialCoordinate(mesh)

    W = FunctionSpace(mesh, "CG", 1)
    if BC is "Dirichlet":
        bcs = [DirichletBC(W, Constant(1.0), 4), DirichletBC(W, Constant(0.0), 3)] 
    else:
        bcs = None
    
    u = Function(W)
    v = TestFunction(W)

    k_eff = as_tensor(k_eff)
    a = inner(k_eff * grad(u), grad(v)) * dx 
    L =  f(x[0]) * v * dx
    if BC is "Neumann":
        a = a + u * v * dx
        bcs = None

    solve(a - L == 0, u, bcs, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
    
    return (u)

