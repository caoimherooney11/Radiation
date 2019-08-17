from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from full_pythonToGeo import pythonToGeo
from math import pi
from makeFullMesh import makeMesh

# domain parameters
mesh_name = "holes_mesh"
domain_dimensions = [0.2, 0.2]# 1.0]
dim = len(domain_dimensions)
radius = 0.2  # O(delta)
scale = 0.5
delta = 0.1
holes_x = int(domain_dimensions[0]/delta)
holes_y = int(domain_dimensions[1]/delta)
num = holes_x * holes_y

#xx = 1 # exponent of nonlinearity ie u**xx
Dirichlet = True
num_constraints = num
norms = []
for i in range(3):
    out = File("Output/u_%d.pvd" % i)
    scale = 2**(-i-1)
    mesh = makeMesh(mesh_name, domain_dimensions, delta, radius, scale)
    
    Z = FunctionSpace(mesh, "CG", 1)
    for i in range(num):
        Z *= FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z, name="solution")
    lam = split(z)
    u = lam[0]
    mu = TestFunctions(Z)
    v = mu[0]
    u_ = z.split()[0]
    
    top_label = num+1
    bottom_label = num+2
    sides_label = num+3
    
    x, y = SpatialCoordinate(mesh)
    n = FacetNormal(mesh)
    radius = radius*delta
    vf = 1. / (4 * pi * radius**2)
    k = 1. # problems arise for larger k. 
    c = 1/delta
    tau = delta**(1/3)
    
    uex = y**2 + x**2
    out.write(u_.interpolate(uex))
    u_.assign(0) 
    #u_.interpolate(uex)
    f = uex - div(k * grad(uex))
    g = -inner(grad(uex), n)
    if Dirichlet:
        bcs = [DirichletBC(V, uex, top_label), DirichletBC(V, uex, bottom_label)]
        flux_bdys = ds(top_label) + ds(bottom_label) + ds(sides_label)
    else:
        bcs = None
        flux_bdys = ds(top_label) + ds(bottom_label) + ds(sides_label)
    
    for eps in [0.25, 0.5, 0.75, 1.0]:
        warning("eps = %f" % eps)
        xx = Constant(eps * 4)
        F = u * v * dx + k * inner(grad(u), grad(v)) * dx - f * v * dx + g * v * flux_bdys
        for i in range(num):
            g1 = - k * inner(grad(uex), n) - c * (uex+tau)**xx + c * vf * assemble((uex+tau)**xx*ds(i+1))
            area = assemble(Constant(1) * ds(i+1, domain=mesh))

            F = F + (c * abs(u+tau)**xx - lam[i+1] + g1) * v * ds(i+1)\
                    + (lam[i+1]/area - c * vf * abs(u+tau)**xx) * mu[i+1] * ds(i+1) 

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
    
    
        solve(F == 0, z, bcs=bcs, solver_parameters=sp)
    out.write(u_)
        
    norm_ = norm(uex-u_)
    norms.append(norm_)
    #print(norm(uex-u_))
    
ratios = []
for i in range(len(norms) - 1):
    ratios.append(norms[i]/norms[i+1])
print(norms)
print(ratios)
