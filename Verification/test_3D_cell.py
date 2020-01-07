from firedrake import *
import numpy as np
from makCellMesh_new import makeCellMesh 

path = "HomogOutput/3D/"
mesh_name = "3D_mesh"
cell_scale = 0.2
global_scale = 0.5
radius = 0.1
k = 1.
c = 1. # c in place of lambda
xi = 0.7; L=8.; stb = 5.67e-8
Tmin = 300
dim = 3
nonlinear = True

def vf(u):
    return 1/((2 * u)**(dim-1) * pi)

mesh_size = []
h = []
times = []
for i in range(3):
    scale = 0.1
    global_scale = 2**(-i)
    mesh = makeCellMesh(mesh_name, radius, cell_scale, global_scale, dim, True)
    mesh_size.append(mesh.num_cells())
    h.append(mesh.num_cells()**(-1/3))

    #V = VectorFunctionSpace(mesh, "CG", 1) 
    V = FunctionSpace(mesh, "CG", 1) 
    u = Function(V)
    v = TestFunction(V)
    n = FacetNormal(mesh)
    vf = 1. / ((2 * radius)**(dim-1) * pi)
    k = 5.
    c = 4.
    T = 1.
    tau = 1.
    
    X  = SpatialCoordinate(mesh) 
    if mms:
        out = File("Output/mms_new.pvd")
        #uex = as_vector([cos(2 * pi * X[0]), cos(2 * pi * X[1]), cos(2 * pi * X[2])])
        #uex = as_vector([cos(2 * pi * X[0]), Constant(1.0), Constant(1.0)])
        uex = cos(2 * pi * X[0])
        out.write(u.interpolate(uex))
        u.assign(0)
        f = uex - div(grad(uex))
        #g = - dot(grad(uex) + Identity(dim), n) - (4 * c * (T + tau)**3 / k) * (uex + X \
        #        - vf * as_vector([assemble((uex[0] + X[0]) * ds(1)), assemble((uex[1] + X[1]) * ds(1)), assemble((uex[2] + X[2]) * ds(1))]))
        g = dot(grad(uex), n)

    #F = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) * dx \
    #        + (inner(grad(u[0]), grad(v[0])) + inner(grad(u[1]), grad(v[1])) \
    #            + inner(grad(u[2]), grad(v[2]))) * dx \
    #        - inner(g, v) * ds(1) - inner(f, v) * dx 
    F = u * v * dx + inner(grad(u), grad(v)) * dx - g * v * ds(1) - f * v * dx
    if not periodic:
        F = F - g * v * ds(2)
    
    sp = {
        "snes_monitor": None,
        #"ksp_monitor": None,
        "mat_type": "matfree",
        "ksp_type": "gmres",
        }
    
    solve(F == 0, u, bcs=None, solver_parameters=sp)
    out.write(u)
    
    if mms:
        norm_ = norm(uex-u)
        norms.append(norm_)
        print(norm_)

    t2 = time.time()
    times.append(t2-t1)

ratios = []
mesh_ratios = []
convergence = []
for i in range(len(norms) - 1):
    ratios.append(norms[i+1] / norms[i])
    mesh_ratios.append(mesh_size[i+1] / mesh_size[i])
    convergence.append(np.log(norms[i+1]/norms[i]) / np.log(h[i+1]/h[i]))

print("times: ", times)
print("mesh size: ", mesh_size)
print("mesh size ratios: ", mesh_ratios)
print("error: ", norms)
print("error ratios: ", ratios)
print("convergence rate: ", convergence)





