from firedrake import Mesh 
from pythonToGeo import pythonToGeo
from subprocess import call 

def makeMesh(mesh_name, domain_dimensions, radius, scale):
    dim = len(domain_dimensions)
    a = pythonToGeo(mesh_name, domain_dimensions, radius)
    call(["gmsh", "-%i" % dim, "%s.geo" % mesh_name, "-clscale", "%f" % scale, "-format", "msh2", "-2"]);

    return Mesh("%s.msh" % mesh_name, dim=dim)
