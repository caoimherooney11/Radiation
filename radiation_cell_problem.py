from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from cell_pythonToGeo import pythonToGeo
from math import pi
from makeCellMesh import makeCellMesh

def solveCellProblem(mesh, radius, k, T, tau, c, VF, nonlinear):
    Z = VectorFunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) * FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z, name = "solution")
    u, lam0, lam1 = split(z)
    v, mu0, mu1 = TestFunctions(Z)
    u_ = z.split()[0]
    n = FacetNormal(mesh)
    #vf = 1. / (4 * pi * radius**2)
    vf = VF(radius)
    
    X  = SpatialCoordinate(mesh) 
    
    area = assemble(Constant(1) * ds(1, domain=mesh))
    if nonlinear:
        F = (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1]))) * dx \
                + (4 * c * (T + tau)**3 / k) * (u[0] + X[0] - lam0) * v[0] * ds(1) \
                + (4 * c * (T + tau)**3 / k) * (u[1] + X[1] - lam1) * v[1] * ds(1) \
                + (lam0/area - vf * (u[0] + X[0])) * mu0 * ds(1) \
                + (lam1/area - vf * (u[1] + X[1])) * mu1 * ds(1) \
                + inner(v, n) * ds(1) 
    else:
        F = (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1]))) * dx \
                + (c / k) * (u[0] + X[0] - lam0) * v[0] * ds(1) \
                + (c / k) * (u[1] + X[1] - lam1) * v[1] * ds(1) \
                + (lam0/area - vf * (u[0] + X[0])) * mu0 * ds(1) \
                + (lam1/area - vf * (u[1] + X[1])) * mu1 * ds(1) \
                + inner(v, n) * ds(1) 
    
    num_constraints = 2
    sp = {
        "snes_monitor": None,
        #"ksp_monitor": None,
        "mat_type": "matfree",
        "ksp_type": "gmres",
        "pc_type": "fieldsplit",
        "pc_fieldsplit_type": "schur",
        "pc_fieldsplit_schur_factorization_type": "diag",
        "pc_fieldsplit_0_fields": "0",
        "pc_fieldsplit_1_fields": ",".join(["%i" % (i+1) for i in range(num_constraints)]),
        "fieldsplit_0": {"ksp_type": "preonly","pc_type": "python", "pc_python_type": 
                         "firedrake.AssembledPC", "assembled_pc_type": "lu", 
                         "assembled_pc_factor_mat_solver_type": "mumps"},
        "fieldsplit_1": {"ksp_type": "preonly","pc_type": "none",},
        }
    
    
    solve(F == 0, z, bcs=None, solver_parameters=sp)
    return u_
    
