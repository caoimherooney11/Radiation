from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
from full_pythonToGeo import pythonToGeo
from math import pi
from makeFullMesh import makeMesh
import time

def solve_direct(mesh_name, domain_dimensions, k, delta, radius, tau, c, scale, BC, f, nonlinear, MMS):
    dim = len(domain_dimensions)
    holes_x = int(domain_dimensions[0]/delta)
    holes_y = int(domain_dimensions[1]/delta)
    num = holes_x * holes_y
    num_constraints = num

    t0 = time.time()
    warning("making direct mesh...")
    mesh = makeMesh(mesh_name, domain_dimensions, delta, radius, scale)
    t1 = time.time()
    warning("direct mesh made")
    mesh_size = mesh.num_cells()
    particle_mesh = mesh_size/num
    mesh_time = t1-t0
    x = SpatialCoordinate(mesh)
    
    warning("solving direct problem...")
    t2 = time.time()
    Z = FunctionSpace(mesh, "CG", 2)
    for i in range(num):
        Z *= FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z, name="direct")
    lam = split(z)
    u = lam[0]
    mu = TestFunctions(Z)
    v = mu[0]
    u_ = z.split()[0]

    top_label = num+1
    bottom_label = num+2
    sides_label = num+3
    
    n = FacetNormal(mesh)
    radius = radius*delta
    vf = 1. / (4 * pi * radius**2)

    if MMS:
        #uex = x[0]**2 + x[1]**2
        uex = x[1]
        f = -div(k * grad(uex))
        g = -inner(grad(uex), n)
        if BC is "Dirichlet":
            bcs = [DirichletBC(V, uex, top_label), DirichletBC(V, uex, bottom_label)] 
            flux_bdys = ds(sides_label)
        else:
            bcs = None
            flux_bdys = ds(top_label) + ds(bottom_label) + ds(sides_label)
    else:
        g = 1e-10
        if BC is "Dirichlet":
            bcs = [DirichletBC(V, Constant(1.0), top_label), DirichletBC(V, Constant(0.0), bottom_label)] 
            flux_bdys = ds(sides_label)
        else:
            bcs = None
            flux_bdys = ds(top_label) + ds(bottom_label) + ds(sides_label)

    u_.assign(0) 
    if nonlinear:
        if MMS:
            N = 11
        else:
            N = 1000
        eps_list = np.linspace(0.0, 1.0, N)

        for eps in eps_list:
            warning("solving for eps = %f" % eps)
            xx = Constant(eps * 4)
            F = k * inner(grad(u), grad(v)) * dx - f(x[0]) * v * dx + g * v * flux_bdys
            for i in range(num):
                area = assemble(Constant(1) * ds(i+1, domain=mesh))
                #warning("area = %f" % area)
                if MMS:
                    g1 = - k * inner(grad(uex), n) - (c/delta) * abs(uex + tau)**xx + (c/delta) * vf * assemble(abs(uex + tau)**xx * ds(i+1))
                else:
                    g1 = Constant(1e-10)
                F = F + ( (c/delta) * abs(u + tau)**xx - lam[i+1] + g1) * v * ds(i+1)\
                        + ( lam[i+1]/area - (c/delta) * vf * abs(u + tau)**xx ) * mu[i+1] * ds(i+1) 
        
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

    else:
        F = k * inner(grad(u), grad(v)) * dx - f(x[0]) * v * dx + g * v * flux_bdys
        for i in range(num):
            area = assemble(Constant(1) * ds(i+1, domain=mesh))
            if MMS:
                g1 = - k * inner(grad(uex), n) - (c/delta) * uex + (c/delta) * vf * assemble(uex * ds(i+1))
            else:
                g1 = Constant(1e-10)
            F = F + ( (c/delta) * u - lam[i+1] + g1) * v * ds(i+1, domain=mesh)\
                    + ( lam[i+1]/area - (c/delta) * vf * u ) * mu[i+1] * ds(i+1, domain=mesh) 
    
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

    if MMS:
        norm_ = norm(uex-u_)
        #norms.append(norm_)
        print(norm(uex-u_))
    else:
        norm_ = 0.0
    
    t3 = time.time() 
    warning("direct problem solved")
    soln_time = t3-t2

    return (u_, mesh_size, particle_mesh, mesh_time, soln_time, norm_)
        
