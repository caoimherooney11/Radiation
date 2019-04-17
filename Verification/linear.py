from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from pythonToGeo import pythonToGeo
from math import pi
from makeMesh import makeMesh

# domain parameters
mesh_name = "single_hole"
domain_dimensions = [1.0, 2.0]# 1.0]
dim = len(domain_dimensions)
radius = 0.2  # size of void within
scale = 0.5
holes_x = 1
holes_y = 2
num = holes_x * holes_y

withIntegral = True
mms = True
num_constraints = num
norms = []
for i in range(4):
    out = File("Output/u.pvd")
    scale = 2**(1-i)
    mesh = makeMesh(mesh_name, domain_dimensions, radius, holes_x, holes_y, scale)
    
    if withIntegral:
        if num_constraints is 1:
            Z = FunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) 
            V = Z.sub(0)
            z = Function(Z)
            u, lam = split(z)
            v, mu = TestFunctions(Z)
        elif num_constraints is 2:
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
        f = uex - div(k * grad(uex))
        g = -inner(grad(uex), n)
    
        if withIntegral:
            if num_constraints is 1:
                g1 = - k * inner(grad(uex), n) - sigma * uex + sigma * vf * assemble(uex*ds(1))
            elif num_constraints is 2:
                g1 = - k * inner(grad(uex), n) - sigma * uex + sigma * vf * assemble(uex*ds(1))
                g2 = - k * inner(grad(uex), n) - sigma * uex + sigma * vf * assemble(uex*ds(2))
        else:
            g1 = - k * inner(grad(uex), n) - sigma * uex
    else:
        f = Constant(0, domain=mesh) 
        g = Constant(0, domain=mesh)
    
    if withIntegral:
        area1 = assemble(Constant(1) * ds(1, domain=mesh)) 
        if num_constraints is 1:
            F = u * v * dx + k * inner(grad(u), grad(v)) * dx - f * v * dx + g * v * flux_bdys\
               + (sigma * u - lam + g1) * v * ds(1) \
               + (lam/area1 - sigma * vf * u) * mu * ds(1) 
        elif num_constraints is 2:
            area2 = assemble(Constant(1) * ds(2, domain=mesh)) 
            F = u * v * dx + k * inner(grad(u), grad(v)) * dx - f * v * dx + g * v * flux_bdys\
               + (sigma * u - lam1 + g1) * v * ds(1) + (sigma * u - lam2 + g2) * v * ds(2)\
               + (lam1/area1 - sigma * vf * u) * mu1 * ds(1) + (lam2/area2 - sigma * vf * u) * mu2 * ds(2) 

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
    else:
        if num_constraints is 1:
            inner_bdys = ds(1)
        elif num_constraints is 2:
            inner_bdys = ds(1) + ds(2)
    
        F = u * v * dx + k * inner(grad(u), grad(v)) * dx - f * v * dx + g * v * flux_bdys\
                + (sigma * u + g1) * v * inner_bdys
        sp = {
                "snes_monitor": None,
                #"ksp_monitor": None,
                "mat_type": "aij",
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps"
                }
    
    
    solve(F == 0, z, bcs=None, solver_parameters=sp)
    out.write(u_)
    
    if mms:
        norm_ = norm(uex-u_)
        norms.append(norm_)
        #print(norm(uex-u_))
    
print(norms)
