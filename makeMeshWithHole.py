from firedrake import *
from cell_pythonToGeo import pythonToGeo
from torusToSquare import torusToSquare
from subprocess import call

mesh_name = "cell_mesh"
hole_scale = 0.5
global_scale = 1.0
pythonToGeo(mesh_name, hole_scale)
call(["gmsh", "-%i" % 2, "%s.geo" % mesh_name, "-clscale", "%f" % global_scale, "-format", "msh2"]);
mesh = Mesh("%s.msh" % mesh_name, dim=3)


#mesh = Mesh("untitled.msh", dim=3)
V = FunctionSpace(mesh, "CG", 1)
v = Function(V)
v.interpolate(Constant(1.0))
File("testing.pvd").write(v)

mesh_new = torusToSquare(mesh)
V = FunctionSpace(mesh_new, "CG", 1)
v = Function(V)
v.interpolate(Constant(1.0))
File("testing.pvd").write(v)

T = mesh_new.coordinates
V = T.function_space()
V = VectorFunctionSpace(mesh_new, "CG", 1)
X = SpatialCoordinate(mesh_new)
u = Function(V)
v = TestFunction(V)
F = inner(grad(u), grad(v)) * dx + inner(u, v) * dx + inner(Constant((0, 0)), v) * dx
bcs = [DirichletBC(V, X-0.1 * X/sqrt(inner(X,X)), 1, "geometric"), 
       DirichletBC(V, 0, [2, 3], "geometric")]
solve(F==0, u, bcs=bcs)
File("u.pvd").write(u)


mesh = Mesh(T)
mesh.coordinates.interpolate(mesh.coordinates - u)
File("round.pvd").write(mesh.coordinates)


