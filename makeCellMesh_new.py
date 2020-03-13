from firedrake import *
from radiation_pythonToGeo_new import pythonToGeo
from subprocess import call
from glueplex_new import make_periodic_mesh

def makeCellMesh(mesh_name, radius, cell_scale, global_scale, dim, periodic):

    pythonToGeo(mesh_name, radius, cell_scale, dim)
    #global_scale = 0.5
    call(["gmsh", "-%i" % dim, "%s.geo" % mesh_name, "-clscale", "%f" % global_scale, "-format", "msh2"]);
    mesh = Mesh("%s.msh" % mesh_name, dim=dim, reorder=False)

    if periodic:
        class Mapping():
        
            def is_slave(self, x):
                eps = 1e-8
                if dim is 2:
                    return (x[0] > 1-eps or x[1] > 1-eps)
                elif dim is 3:
                    return (x[0] > 1-eps or x[1] > 1-eps or x[2] > 1-eps)

            def is_master(self, x):
                eps = 1e-8
                if dim is 2:
                    return (x[0] < eps and x[1] <= 1-eps) or \
                        (x[1] < eps and x[0] <= 1-eps) 
                if dim is 3:
                    return (x[0] < eps and x[1] <= 1-eps and x[2] <= 1-eps) or \
                        (x[1] < eps and x[0] <= 1-eps and x[2] <= 1-eps) or \
                        (x[2] < eps and x[0] <= 1-eps and x[1] <= 1-eps) 
        
            def map_to_master(self, x):
                master = x.copy()
                eps = 1e-8
                if x[0] > 1-eps:
                    master[0] -= 1
                if x[1] > 1-eps:
                    master[1] -= 1
                if dim is 3:
                    if x[2] > 1-eps:
                        master[2] -= 1
                return master
        
        mapping = Mapping()
        periodic_mesh = make_periodic_mesh(mesh, mapping)

    else: 
        periodic_mesh = mesh
    
    W = VectorFunctionSpace(periodic_mesh, "CG", 1)
    w = Function(W).interpolate(SpatialCoordinate(periodic_mesh))
    File("test.pvd").write(w)

    #mesh_new = torusToSquare(mesh)
    #T = mesh_new.coordinates
    #V = T.function_space()
    #V = VectorFunctionSpace(mesh_new, "CG", 1)
    #X = SpatialCoordinate(mesh_new)
    #u = Function(V)
    #v = TestFunction(V)
    #F = inner(grad(u), grad(v)) * dx + inner(u, v) * dx + inner(Constant((0, 0)), v) * dx
    #bcs = [DirichletBC(V, X - radius * X/sqrt(inner(X,X)), 1, "geometric"), 
    #       DirichletBC(V, 0, [2, 3], "geometric")]
    #solve(F==0, u, bcs=bcs)
    ##File("u.pvd").write(u)

    #mesh = Mesh(T)
    #mesh.coordinates.interpolate(mesh.coordinates - u)
    #File("round.pvd").write(mesh.coordinates)

    return periodic_mesh


