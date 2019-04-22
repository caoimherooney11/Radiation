from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from full_pythonToGeo import pythonToGeo
from math import pi
from makeFullMesh import makeMesh

# domain parameters
mesh_name = "holes_mesh"
domain_dimensions = [2.0, 2.0]# 1.0]
dim = len(domain_dimensions)
radius = 0.2  # size of void within
scale = 0.5
holes_x = 2
holes_y = 2
num = holes_x * holes_y

#xx = 1 # exponent of nonlinearity ie u**xx
mms = True
num_constraints = num
norms = []
for i in range(4):
    out = File("Output/u.pvd")
    scale = 2**(1-i)
    mesh = makeMesh(mesh_name, domain_dimensions, radius, holes_x, holes_y, scale)
    
    Z = FunctionSpace(mesh, "CG", 1)
    for i in range(num):
        Z *= FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z)
    lam = split(z)
    u = lam[0]
    mu = TestFunctions(Z)
    v = mu[0]
    u_ = z.split()[0]
    
    top_label = num+1
    bottom_label = num+2
    sides_label = num+3
    flux_bdys = ds(top_label) + ds(bottom_label) + ds(sides_label)
    
    x, y = SpatialCoordinate(mesh)
    n = FacetNormal(mesh)
    vf = 1. #/ (4 * pi * radius**2)
    k = 1.
    sigma = 1.
    
    if mms:
        uex = y**2 + x**2
        out.write(u_.interpolate(uex))
        u_.assign(0) 
        #u_.interpolate(uex)
        f = uex - div(k * grad(uex))
        g = -inner(grad(uex), n)
    
    else:
        f = Constant(0, domain=mesh) 
        g = Constant(0, domain=mesh)
    
    for eps in [0.25, 0.5, 0.75, 1.0]:
        xx = eps * 4
        F = u * v * dx + k * inner(grad(u), grad(v)) * dx - f * v * dx + g * v * flux_bdys
        for i in range(num):
            g1 = - k * inner(grad(uex), n) - sigma * uex**xx + sigma * vf * assemble(uex**xx*ds(i+1))
            area = assemble(Constant(1) * ds(i+1, domain=mesh))

            F = F + (sigma * u**xx - lam[i+1] + g1) * v * ds(i+1)\
                    + (lam[i+1]/area - sigma * vf * u**xx) * mu[i+1] * ds(i+1) 

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
            "fieldsplit_0": {"ksp_type": "preonly","pc_type": "python", "pc_python_type": "firedrake.AssembledPC", "assembled_pc_type": "lu", "assembled_pc_factor_mat_solver_type": "mumps"},
            "fieldsplit_1": {"ksp_type": "preonly","pc_type": "none",},
            }
    
    
        solve(F == 0, z, bcs=None, solver_parameters=sp)
    out.write(u_)
        
    if mms:
        norm_ = norm(uex-u_)
        norms.append(norm_)
        #print(norm(uex-u_))
    
print(norms)
