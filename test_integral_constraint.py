from firedrake import *

for i in range(5):
    print(2**i)
    out = File("Output/test.pvd")
    mesh = UnitSquareMesh(2**i, 2**i)
    Z = FunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "R", 0) #* FunctionSpace(mesh, "R", 0)
    V = Z.sub(0)
    z = Function(Z)
    u, lam = split(z)
    #u, lam, lam_ = split(z)
    v, mu  = TestFunctions(Z)
    #v, mu, mu_  = TestFunctions(Z)
    u_ = z.split()[0]
    
    bcs = [DirichletBC(V, Constant('1.0'), 4), DirichletBC(V, Constant('0.0'), 3)]
    
    x, y = SpatialCoordinate(mesh)
    uex = y**2
    out.write(u_.interpolate(uex))
    n = FacetNormal(mesh)
    f = -div(grad(uex))
    g2 = assemble(uex * ds(2))
    g1 = - inner(grad(uex), n)
    #g1 = assemble(uex * ds(1))
    
    area1 = assemble(Constant(1) * ds(1, domain=mesh))
    area2 = assemble(Constant(1) * ds(2, domain=mesh))
    F = inner(grad(u), grad(v)) * dx - f * v * dx \
            - inner(grad(u), n) * v * ds(2) + (lam/area2 - u + g2) * mu * ds(2) \
            + g1 * v * ds(1)
            #- inner(grad(u), n) * v * ds(1) + (lam_/area1 - u + g1) * mu_ * ds(1)
    sp = {
        "snes_monitor": None,
        "ksp_monitor": None,
    }
    solve(F == 0, z, bcs=bcs, solver_parameters=sp)
    
    out.write(u_)
    print(norm(uex-u_))
