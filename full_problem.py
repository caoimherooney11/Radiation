from firedrake import *
import numpy as np
from interpolate_keff import keff

def solve_full(domain_dimensions, lists):
    nx = domain_dimensions[0] 
    ny = domain_dimensions[1]
    n = 10
    mesh = RectangleMesh(10 * nx, 10 * ny, nx, ny)
    
    P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, P1)
    V = TensorFunctionSpace(mesh, P1)
    bcs = [DirichletBC(W, Constant(1.0), 4), DirichletBC(W, Constant(0.0), 3)] 
    
    u = Function(W)
    v = TestFunction(W)
    du = Function(W)
    g = Constant(0.0)
    k_eff = Function(V)
    dk_eff = Function(V)
    
    F = inner(k_eff * grad(u), grad(v)) * dx
    J = inner(k_eff * grad(du), grad(v)) * dx \
            + inner(du * dk_eff * grad(u), grad(v)) * dx 
    
    uvec = u.vector()
    max_iter = 3
    for i in range(max_iter):
        temp1 = k_eff.vector().get_local()
        temp2 = dk_eff.vector().get_local()
        tempu = u.vector().get_local()
        for j in range(len(tempu)):
            for i in range(4):
                temp1[4*j + i] = keff(tempu[j], lists, entry=i+1)
                temp2[4*j + i] = keff(tempu[j], lists, entry=i+1, derivative=1)
        k_eff.vector().set_local(temp1)
        dk_eff.vector().set_local(temp2)
    
        solve(J + F == 0, du, bcs, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})
        bcs = [DirichletBC(W, Constant(0.0), 4), DirichletBC(W, Constant(0.0), 3)] 
        uvec += du.vector()
    
    File("Output/T.pvd").write(u)
    
    return (u)

