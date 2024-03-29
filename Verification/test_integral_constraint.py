from firedrake import *
from makeMesh import makeMesh

num_constraints = 1
holes = True
for j in range(1):
    out = File("Output/test.pvd")
    i = j + 5
    print(2**i)
    if holes:
        if num_constraints is 1:
            #mesh = Mesh("single_hole.msh", dim=2)
            mesh = makeMesh()
        else:
            mesh = Mesh("many_holes.msh", dim=2)
    else:
        mesh = UnitSquareMesh(2**i, 2**i)

    if num_constraints is 1:
        Z = FunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0)
        V = Z.sub(0)
        z = Function(Z)
        u, lam = split(z)
        v, mu  = TestFunctions(Z)
    if num_constraints is 2:
        Z = FunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) * FunctionSpace(mesh, "R", 0)
        V = Z.sub(0)
        z = Function(Z)
        u, lam, lam_ = split(z)
        v, mu, mu_  = TestFunctions(Z)
    u_ = z.split()[0]
    u_.rename("output")
    
    #bcs = [DirichletBC(V, Constant('1.0'), 4), DirichletBC(V, Constant('0.0'), 3)]
    
    x, y = SpatialCoordinate(mesh)
    uex = y
    out.write(u_.interpolate(uex))
    u_.assign(0)
    n = FacetNormal(mesh)
    f = uex - div(grad(uex))
    g = -inner(grad(uex), n)
    g1 = 2*assemble(uex * ds(1))
    print(g1)
    area1 = assemble(Constant(1) * ds(1, domain=mesh))
    print(area1)

    if num_constraints is 1:
        flux_bdys = ds(2) + ds(3) + ds(4)
        F = u * v * dx + inner(grad(u), grad(v)) * dx - f * v * dx + (g * v) * flux_bdys\
            -inner(grad(u), n) * v * ds(1)# + (lam/area1 - u + g1/area1) * mu * ds(1)
        #F = u * v * dx + inner(grad(u), grad(v)) * dx - lam * v * ds(1) - mu * (u-g/area1)*ds(1) \
        #    - f * v * dx + (g * v) * flux_bdys - inner(grad(u), n) * v * ds(1)
    if num_constraints is 2:
        if holes:
            flux_bdys = ds(3) + ds(4) + ds(5)
            g2 = assemble(uex * ds(2))
            area2 = assemble(Constant(1) * ds(2, domain=mesh))
            F = u * v * dx + inner(grad(u), grad(v)) * dx - f * v * dx + (g * v) * flux_bdys\
                + (lam/area1 - u + g1/area1) * mu * ds(1)\
                + (lam_/area2 - u + g2/area2) * mu_ * ds(2)
        else:
            flux_bdys = ds(2) + ds(4)
            g3 = assemble(uex * ds(3))
            area3 = assemble(Constant(1) * ds(3, domain=mesh))
            F = u * v * dx + inner(grad(u), grad(v)) * dx - f * v * dx + (g * v) * flux_bdys\
                + (lam/area1 - u + g1/area1) * mu * ds(1)\
                + (lam_/area3 - u + g3/area3) * mu_ * ds(3)

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
    print(assemble(u_ * ds(1)))
    
    out.write(u_)
    print("||u-uex||", norm(uex-u_))
