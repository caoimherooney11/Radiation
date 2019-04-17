from firedrake import Mesh 
from full_pythonToGeo import pythonToGeo
from subprocess import call 

def makeMesh(mesh_name, domain_dimensions, radius, holes_x, holes_y, scale):
    dim = len(domain_dimensions)
    pythonToGeo(mesh_name, domain_dimensions, radius, holes_x, holes_y)
    call(["gmsh", "-%i" % dim, "%s.geo" % mesh_name, "-clscale", "%f" % scale, "-format", "msh2", "-2"]);

    return Mesh("%s.msh" % mesh_name, dim=dim)
 
