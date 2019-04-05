from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from pythonToGeo import pythonToGeo
from math import pi
from makeMesh import makeMesh

# domain parameters
(mesh_name,domain_dimensions, radius, num) = makeMesh()
#domain_dimensions = [1.0, 1.0]# 1.0]
dim = len(domain_dimensions)
#radius = 0.2  # size of void within
#scale = 1.0
#pythonToGeo(mesh_name, domain_dimensions, radius)
mesh = Mesh("%s.msh" % mesh_name, dim=dim)

withIntegral = True
mms = True
if withIntegral:
    Z = FunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) * FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z)
    u, lam1, lam2 = split(z)
    v, mu1, mu2 = TestFunctions(Z)
    u_ = z.split()[0]
else:
    V = FunctionSpace(mesh, "CG", 1)
    z = Function(V)
    u = z
    v = TestFunction(V)
    u_ = z

# dirichlet boundary labels: top  = number of #holes+1, bottom = #holes+2
top_label = num+1
bottom_label = num+2
bcs = [DirichletBC(V, Constant('2.0'), top_label), DirichletBC(V, Constant('0.0'), bottom_label)]

sp = {
    "snes_monitor": None,
    # "ksp_type": "gmres",
    # "pc_type": "hypre",
    "ksp_monitor": None
}


if dim == 2:
    x, y = SpatialCoordinate(mesh)
else: 
    x, z_, y = SpatialCoordinate(mesh) # rotating axis so that mms works in both 2D and 3D
n = FacetNormal(mesh)
vf = 1 / (4 * pi * radius**2)
k = 1.
sigma = 1.

#out = File("/home/wechsung/Dropbox/temp/u.pvd")
out = File("Output/u.pvd")
if mms:
    uex = y
    f = -div(k * grad(uex))
    eps = Constant(1)
    if withIntegral:
        #g = - k * inner(grad(uex), n) - sigma * uex**4# + sigma * vf * assemble(eps * uex**4*ds(4))
        #for i in range(num):
        #    g = g + sigma * vf * assemble(eps * uex**4 * ds(i+1))
        g1 = - k * inner(grad(uex), n) - sigma * uex**4 + sigma * vf * assemble(eps * uex**4*ds(1))
        g2 = - k * inner(grad(uex), n) - sigma * uex**4 + sigma * vf * assemble(eps * uex**4*ds(2))
    else:
        g = - k * inner(grad(uex), n) - sigma * uex**4 
    out.write(u_.interpolate(uex))
    u_.assign(0) 
else:
    f = Constant(0, domain=mesh) 
    g = Constant(0, domain=mesh)

if withIntegral:
   # F = k * inner(grad(u), grad(v)) * dx - f*v*dx #+ (sigma * u**4 - lam + g) * v * ds(4) + (lam/area - sigma * vf * u**4) * mu * ds(4)
   # for i in range(num):
   #     area = assemble(Constant(1) * ds(i+1, domain=mesh)) 
   #     F = F + (sigma * u**4 - lam + g) * v * ds(i+1) + (lam/area - sigma * vf * u**4) * mu * ds(i+1)
   area1 = assemble(Constant(1) * ds(1, domain=mesh)) 
   area2 = assemble(Constant(1) * ds(2, domain=mesh)) 
   F = k * inner(grad(u), grad(v)) * dx - f * v * dx \
           + (sigma * u**4 - lam1 + g1) * v * ds(1) + (sigma * u**4 - lam2 + g2) * v * ds(2)\
           + (lam1/area1 - sigma * vf * u**4) * mu1 * ds(1) + (lam2/area2 - sigma * vf * u**4) * mu2 * ds(2) 
else:
    F = k * inner(grad(u), grad(v)) * dx - f * v * dx #+ (sigma * u**4 + g) * v * ds(4)
    for i in range(num):
        print(i+1)
        F = F + (sigma * u**4 + g) * v * ds(i+1)

solve(F == 0, z, bcs=bcs, solver_parameters=sp)
out.write(u_)

if mms:
    print(norm(uex-u_))
