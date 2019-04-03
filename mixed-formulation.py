from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from pythonToGeo import pythonToGeo
from makeMesh import makeMesh
from math import pi

# domain parameters
mesh_name = "new_test"
domain_dimensions = [1.0, 1.0]  # unit square
radius = 0.2  # size of void within
scale = 0.5
pythonToGeo(mesh_name, domain_dimensions, radius)
mesh = Mesh("%s.msh" % mesh_name, dim=2)

withIntegral = False
mms = True
if withIntegral:
    Z = FunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z)
    u, lam = split(z)
    v, mu = TestFunctions(Z)
    u_ = z.split()[0]
else:
    V = FunctionSpace(mesh, "CG", 1)
    z = Function(V)
    u = z
    v = TestFunction(V)
    u_ = z

bcs = [DirichletBC(V, Constant('1.0'), 1), DirichletBC(V, Constant('0.0'), 2)]

sp = {
    "snes_monitor": None,
    # "ksp_type": "gmres",
    # "pc_type": "hypre",
    "ksp_monitor": None
}

area = assemble(Constant(1) * ds(4, domain=mesh))

x, y = SpatialCoordinate(mesh)
n = FacetNormal(mesh)

out = File("/home/wechsung/Dropbox/temp/u.pvd")
if mms:
    uex = y
    f = -div(grad(uex))
    eps = Constant(1)
    if withIntegral:
        g = uex**4 - inner(grad(uex), n) - assemble(eps * uex**4*ds(4))
    else:
        g = uex**4 - inner(grad(uex), n)
    out.write(u_.interpolate(uex))
    u_.assign(0)
else:
    f = Constant(0, domain=mesh)
    g = Constant(0, domain=mesh)

if withIntegral:
    F = inner(grad(u), grad(v)) * dx - f*v*dx - (u**4 - lam - g)*v*ds(4) - (lam/area - u**4) * mu * ds(4)
else:
    F = inner(grad(u), grad(v)) * dx - f*v*dx - (u**4 - g)*v*ds(4)

solve(F == 0, z, bcs=bcs, solver_parameters=sp)
out.write(u_)

if mms:
    print(norm(uex-u_))
