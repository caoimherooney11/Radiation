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
    a = k * inner(grad(u1), grad(v1)) * dx - eps * \
        sigma * u1**4 * v1 * ds(4) + f * v1 * dx
    sp = {
        "snes_monitor": None,
        "ksp_type": "gmres",
        "pc_type": "hypre",
        # "ksp_monitor": None
    }
    solve(a == 0, u1, bcs=bcs, solver_parameters=sp)
    T_output.write(u1)


# manual Newton
u = Function(V)
for bc in bcs:
    bc.apply(u)
v = TestFunction(V)
du = TrialFunction(V)
uu = Function(V)

F = k * inner(grad(u), grad(v)) * dx - sigma * u**4 * v * ds(4)
J = k * inner(grad(du), grad(v)) * dx - 4 * sigma * u**3 * du * v * ds(4)

uvec = u.vector()
max_iter = 10
error_list = []
max_tol = 1e-8
error = 1.0
eps = Constant(1e-1)
stepsize = 1.0
bcs = [DirichletBC(V, Constant(0.0), 1), DirichletBC(V, Constant(0.0), 2)]
# for i in range(max_iter):
while error > max_tol:
    tempu = u.vector().get_local()
    A = assemble(J)
    b = assemble(-F)
    # A.force_evaluation()
    # import IPython; IPython.embed()

    temp = Function(V)
    tempvec = temp.vector()
    # build vector containing int_gamma 4*u**3 test_j ds
    firstintegral = Function(V)
    firstintegralvec = firstintegral.vector()
    for i in range(0, len(tempvec)):
        tempvec[i] = 1.
        firstintegralvec[i] = assemble(4 * u **3 * temp * ds(4))
        tempvec[i] = 0.

    # build vector containing int_gamma trial_i ds
    secondintegral = Function(V)
    secondintegralvec = secondintegral.vector()
    for i in range(0, len(tempvec)):
        tempvec[i] = 1.
        secondintegralvec[i] = assemble(eps * temp * ds(4))
        tempvec[i] = 0.

    # Now modify the matrix
    Ap = A.petscmat
    Ap.setOption(Ap.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    Ap.setOption(Ap.Option.NEW_NONZERO_LOCATION_ERR, False)
    Ap.setOption(Ap.Option.UNUSED_NONZERO_LOCATION_ERR, False)
    for i in range(0, len(tempvec)):
        for j in range(0, len(tempvec)):
            val = -(sigma / (4 * pi * radius**2)) * firstintegralvec[j] * secondintegralvec[i]
            if abs(val) > 1e-10:
                Ap.setValue(i, j, val, addv=PETSc.InsertMode.ADD_VALUES)

    Ap.assemble()

    # Now modify the rhs
    bvec = b.vector()
    for i in range(0, len(tempvec)):
        tempvec[i] = 1.
        bvec[i] = bvec[i] + (sigma/(4*pi*radius**2)) * assemble(u ** 4 * ds(4)) * assemble(eps * temp * ds(4))
        tempvec[i] = 0.

    solve(A, uu, b, bcs=bcs, solver_parameters=sp)

    error = abs(assemble((inner(uu, uu) * dx)) ** 0.5 /
                assemble(inner(u+uu, u+uu) * dx) ** 0.5)

    uvec += stepsize * uu.vector()
    T_manual_output.write(u)

    print(error)
    # error_list.append(error)

# print(error_list)
