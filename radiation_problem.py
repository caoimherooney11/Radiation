from firedrake import *
import numpy as np
from pythonToGeo import pythonToGeo
from makeMesh import makeMesh

# domain parameters
mesh_name = "new_test"
domain_dimensions = [1.0, 1.0] # unit square
radius = 0.2 # size of void within
scale = 0.5
pythonToGeo(mesh_name, domain_dimensions, radius)
#mesh = makeMesh(mesh_name, domain_dimensions, radius, scale)
mesh = Mesh("%s.msh" % mesh_name, dim=2)
#mesh = UnitSquareMesh(10,10)


# linear variational problem
V = FunctionSpace(mesh, "CG", 1)
u1 = Function(V)
v1 = TestFunction(V)
f = Constant(0.0)   
k = 1.0
sigma = 1.0
T_output = File("T.pvd")
T_manual_output = File("T_manual.pvd")

bcs = [DirichletBC(V, Constant('1.0'), 1), DirichletBC(V, Constant('0.0'), 2)]

epsilon = np.linspace(0.01, 0.1, 10)
epsilon = [1.0]
for eps in epsilon:
    print('Solving for eps = ', eps)
    a = k * inner(grad(u1), grad(v1)) * dx - eps * sigma * u1**4 * v1 * ds(4) +  f * v1 * dx
    sp = {"snes_monitor": None, "ksp_type" : "gmres", "pc_type" : "hypre", "ksp_monitor" : True}
    solve(a == 0, u1, bcs=bcs, solver_parameters = sp)
    T_output.write(u1)


# manual Newton
u = Function(V)
v = TestFunction(V)
du = Function(V)

F = k * inner(grad(u), grad(v)) * dx - sigma * u**4 * v * ds(4)
J = k * inner(grad(du), grad(v)) * dx - 4 * sigma * u**3 * du * v * ds(4)

uvec = u.vector()
max_iter = 10 
error_list = []
max_tol = 1e-8
error = 1.0
#for i in range(max_iter):
while error > max_tol:
    tempu = u.vector().get_local()
    solve(J + F == 0, du, bcs=bcs, solver_parameters = {"newton_solver": {"relative_tolerance": 1e-6}})

    error = abs(assemble((inner(du, du) * dx)) ** 0.5 / assemble(inner(u+du, u+du) * dx) ** 0.5)
    print(error)
    bcs = [DirichletBC(V, Constant(0.0), 1), DirichletBC(V, Constant(0.0), 2)]
    uvec += du.vector()
    T_manual_output.write(u)

    #error = abs(assemble((inner(u1 - u, u1 - u) * dx)) ** 0.5 / assemble(inner(u1, u1) * dx) ** 0.5)
    #print(error)
    #error_list.append(error)

#print(error_list)


