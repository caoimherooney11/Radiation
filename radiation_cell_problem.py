from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from cell_pythonToGeo import pythonToGeo
from math import pi
from makeCellMesh import makeCellMesh

def solveCellProblem(mesh, radius, k, T, tau, c, VF, nonlinear, dim):
    Z = VectorFunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) * FunctionSpace(mesh, "R", 0)
    if dim is 3:
        Z = Z * FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z, name = "solution")
    if dim is 2:
        u, lam0, lam1 = split(z)
        v, mu0, mu1 = TestFunctions(Z)
    elif dim is 3:
        u, lam0, lam1, lam2 = split(z)
        v, mu0, mu1, mu2 = TestFunctions(Z)
    u_ = z.split()[0]
    n = FacetNormal(mesh)
    vf = Constant(VF(radius))
    
    X  = SpatialCoordinate(mesh) 
    
    area = assemble(Constant(1) * ds(1, domain=mesh))
    if nonlinear:
        F = (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1]))) * dx \
                + (4 * Constant(c) * Constant(T + tau)**3 ) * (u[0] + X[0] - lam0) * v[0] * ds(1) \
                + (4 * Constant(c) * Constant(T + tau)**3 ) * (u[1] + X[1] - lam1) * v[1] * ds(1) \
                + (lam0/Constant(area) - Constant(vf) * (u[0] + X[0])) * mu0 * ds(1) \
                + (lam1/Constant(area) - Constant(vf) * (u[1] + X[1])) * mu1 * ds(1) \
                + inner(v, n) * ds(1) 
        if dim is 3:
            F = F + inner(grad(u[2]), grad(v[2])) * dx \
                + (4 * Constant(c) * Constant(T + tau)**3 ) * (u[2] + X[2] - lam2) * v[2] * ds(1) \
                + (lam2/Constant(area) - Constant(vf) * (u[2] + X[2])) * mu2 * ds(1) 
    else:
        F = (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1]))) * dx \
                + Constant(c) * (u[0] + X[0] - lam0) * v[0] * ds(1) \
                + Constant(c) * (u[1] + X[1] - lam1) * v[1] * ds(1) \
                + (lam0/area - vf * (u[0] + X[0])) * mu0 * ds(1) \
                + (lam1/area - vf * (u[1] + X[1])) * mu1 * ds(1) \
                + inner(v, n) * ds(1) 
        if dim is 3:
            F = F + inner(grad(u[2]), grad(v[2])) * dx \
                + Constant(c) * (u[2] + X[2] - lam2) * v[2] * ds(2) \
                + (lam2/area - vf * (u[2] + X[2])) * mu2 * ds(2) 
    
    num_constraints = dim
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
    
