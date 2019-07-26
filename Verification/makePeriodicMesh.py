from firedrake import *
from glueplex import make_periodic_mesh
from subprocess import call
from periodic_pythonToGeo import pythonToGeo

def getPeriodicMesh(mesh_name, scale):

    #mesh = UnitSquareMesh(4, 4, reorder=False)
    global_scale = 1.0
    scale1 = 1/scale
    scale2 = scale
    pythonToGeo(mesh_name, 0.25, scale1, scale2)
    call(["gmsh", "-%i" % 2, "%s.geo" % mesh_name, "-clscale", "%f" % global_scale, "-format", "msh2"]);
    mesh = Mesh("%s.msh" % mesh_name, reorder=False)

    class Mapping():
    
        def is_slave(self, x):
            eps = 1e-8
            return (x[0] > 1-eps or x[1] > 1-eps)
    
        def map_to_master(self, x):
            master = x.copy()
            eps = 1e-8
            if x[0] > 1-eps:
                master[0] -= 1
            if x[1] > 1-eps:
                master[1] -= 1
            return master
    
    mapping = Mapping()
    periodic_mesh = make_periodic_mesh(mesh, mapping)
    
    #W = VectorFunctionSpace(periodic_mesh, "CG", 1)
    #w = Function(W).interpolate(SpatialCoordinate(periodic_mesh))
    #File("test.pvd").write(w)
    #
    #import IPython; IPython.embed()

    return periodic_mesh
