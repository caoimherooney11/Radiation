from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from math import pi
from makCellMesh_new import makeCellMesh
import time

mms = True
periodic = True # periodic very slow
dim = 3
norms = []
mesh_size = []
h = []
times = []
for i in range(3):
    t1 = time.time()
    mesh_name = "cell_mesh"
    radius = 0.1
    scale = 0.1
    global_scale = 2**(-i)
    print(scale)
    mesh = makeCellMesh(mesh_name, radius, scale, global_scale, dim, periodic)
    mesh_size.append(mesh.num_cells())
    h.append(mesh.num_cells()**(-1/3))
    
    #V = VectorFunctionSpace(mesh, "CG", 1) 
    #V = FunctionSpace(mesh, "CG", 1) 
    #u = Function(V)
    #v = TestFunction(V)
    Z = VectorFunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) \
            * FunctionSpace(mesh, "R", 0) * FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z)
    u, lam0, lam1, lam2 = split(z)
    v, mu0, mu1, mu2 = TestFunctions(Z)
    u_ = z.split()[0]
    n = FacetNormal(mesh)
    vf = 1. / ((2 * radius)**(dim-1) * pi)
    k = 5.
    c = 4.
    T = 1.
    tau = 1.

    
    X  = SpatialCoordinate(mesh) 
    if mms:
        out = File("Output/mms_new.pvd")
        #uex = as_vector([cos(2 * pi * X[0]), cos(2 * pi * X[1]), cos(2 * pi * X[2])])
        uex = as_vector([cos(2 * pi * X[0]), Constant(1.0), Constant(1.0)])
        out.write(u_.interpolate(uex))
        u_.assign(0)
        f = uex - div(grad(uex))
        g = - dot(grad(uex) + Identity(dim), n) - (4 * c * (T + tau)**3 / k) * (uex + X \
                - vf * as_vector([assemble((uex[0] + X[0]) * ds(1)), assemble((uex[1] + X[1]) * ds(1)), assemble((uex[2] + X[2]) * ds(1))]))
        #g = -dot(grad(uex) + Identity(3), n)

    area = assemble(Constant(1) * ds(1, domain=mesh))
    F = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) * dx \
            + (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1]))+ inner(grad(u[2]), grad(v[2]))) * dx \
            + inner(g, v) * ds(1) + inner(v, n) * ds(1) \
            - inner(f, v) * dx \
            + (4 * c * (T + tau)**3 / k) * (u[0] + X[0] - lam0) * v[0] * ds(1) \
            + (4 * c * (T + tau)**3 / k) * (u[1] + X[1] - lam1) * v[1] * ds(1) \
            + (4 * c * (T + tau)**3 / k) * (u[2] + X[2] - lam2) * v[2] * ds(1) \
            + (lam0/Constant(area) - Constant(vf) * (u[0] + X[0])) * mu0 * ds(1) \
            + (lam1/Constant(area) - Constant(vf) * (u[1] + X[1])) * mu1 * ds(1) \
            + (lam2/Constant(area) - Constant(vf) * (u[2] + X[2])) * mu2 * ds(1) 
    if not periodic:
        F = F - g * v * ds(2)
    
    sp = {
        "snes_monitor": None,
        #"ksp_monitor": None,
        "mat_type": "matfree",
        "ksp_type": "gmres",
        }
    
    solve(F == 0, z, bcs=None, solver_parameters=sp)
    out.write(u_)
    
    if mms:
        norm_ = norm(uex-u_)
        norms.append(norm_)
        print(norm_)

    t2 = time.time()
    times.append(t2-t1)

ratios = []
mesh_ratios = []
convergence = []
for i in range(len(norms) - 1):
    ratios.append(norms[i+1] / norms[i])
    mesh_ratios.append(mesh_size[i+1] / mesh_size[i])
    convergence.append(np.log(norms[i+1]/norms[i]) / np.log(h[i+1]/h[i]))

print("times: ", times)
print("mesh size: ", mesh_size)
print("mesh size ratios: ", mesh_ratios)
print("error: ", norms)
print("error ratios: ", ratios)
print("convergence rate: ", convergence)
    
