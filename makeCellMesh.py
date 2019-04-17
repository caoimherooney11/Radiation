from firedrake import *
from cell_pythonToGeo import pythonToGeo
from torusToSquare import torusToSquare
from subprocess import call

def makeCellMesh(mesh_name, radius, hole_scale, global_scale):
    pythonToGeo(mesh_name, radius, hole_scale)
    call(["gmsh", "-%i" % 2, "%s.geo" % mesh_name, "-clscale", "%f" % global_scale, "-format", "msh2"]);
    mesh = Mesh("%s.msh" % mesh_name, dim=3)

    mesh_new = torusToSquare(mesh)
    T = mesh_new.coordinates
    V = T.function_space()
    V = VectorFunctionSpace(mesh_new, "CG", 1)
    X = SpatialCoordinate(mesh_new)
    u = Function(V)
    v = TestFunction(V)
    F = inner(grad(u), grad(v)) * dx + inner(u, v) * dx + inner(Constant((0, 0)), v) * dx
    bcs = [DirichletBC(V, X - radius * X/sqrt(inner(X,X)), 1, "geometric"), 
           DirichletBC(V, 0, [2, 3], "geometric")]
    solve(F==0, u, bcs=bcs)
    #File("u.pvd").write(u)

    mesh = Mesh(T)
    mesh.coordinates.interpolate(mesh.coordinates - u)
    File("round.pvd").write(mesh.coordinates)

    return mesh


