from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from math import pi
from makeCellMesh import makeCellMesh
from makePeriodicMesh import getPeriodicMesh

mms = True
norms = []
for i in range(4):
    mesh_name = "cell_mesh"
    radius = 0.25
    scale = 2**(-i)/10
    global_scale = 1.0
    print(scale)
    mesh = makeCellMesh(mesh_name, radius, scale, global_scale)
    
    Z = VectorFunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) * FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z)
    u, lam0, lam1 = split(z)
    v, mu0, mu1 = TestFunctions(Z)
    u_ = z.split()[0]
    n = FacetNormal(mesh)
    vf = 1. / (4 * pi * radius**2)
    k = 1.
    sigma = 1.
    T = 1.
    
    X  = SpatialCoordinate(mesh) 
    if mms:
        out = File("Output/mms.pvd")
        uex = as_vector([cos(2 * pi * X[0])*cos(2 * pi * X[1]), sin(2 * pi * X[0])*sin(2 * pi * X[1])])
        out.write(u_.interpolate(uex))
        u_.assign(0)
        f = - k * div(grad(uex))
        g = - k * dot(grad(uex) + Identity(2), n) - 4 * sigma * T**3 * (uex + X \
                - vf * as_vector([assemble((uex[0] + X[0]) * ds(1)), assemble((uex[1] + X[1]) * ds(1))]))
    else:
        out = File("Output/u.pvd")
        f = as_vector([Constant(0.0), Constant(0.0)])
        g = as_vector([Constant(0.0), Constant(0.0)])
    
    
    area = assemble(Constant(1) * ds(1, domain=mesh))
    F = k * (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1]))) * dx - inner(f, v) * dx\
            + (4 * sigma * T / k) * (u[0] + X[0] - lam0) * v[0] * ds(1) \
            + (4 * sigma * T / k) * (u[1] + X[1] - lam1) * v[1] * ds(1) \
            + (lam0/area - vf * (u[0] + X[0])) * mu0 * ds(1) \
            + (lam1/area - vf * (u[1] + X[1])) * mu1 * ds(1) \
            + inner(v, n) * ds(1) \
            + inner(g, v) * ds(1)
    
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
    out.write(u_)
    
    if mms:
        norm_ = norm(uex-u_)
        norms.append(norm_)

print(norms)
    
