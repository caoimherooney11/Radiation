from firedrake import *
import numpy as np
from interpolate_keff import keff

def solve_full(domain_dimensions, k_eff, BC, f):
    nx = domain_dimensions[0] 
    ny = domain_dimensions[1]
    n = 10
    mesh = RectangleMesh(200 * nx, 200 * ny, nx, ny)
    x = SpatialCoordinate(mesh)

    
    P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, P1)
    if BC is "Dirichlet":
        bcs = [DirichletBC(W, Constant(1.0), 4), DirichletBC(W, Constant(0.0), 3)] 
    
    u = Function(W)
    v = TestFunction(W)

    F = inner(k_eff * grad(u), grad(v)) * dx - f * v * dx
    if BC is "Neumann":
        F = F + u * v * dx
        bcs = None

    solve(F == 0, u, bcs, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
    
    return (u)

