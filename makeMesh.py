from firedrake import Mesh 
from pythonToGeo import pythonToGeo
from subprocess import call 

def makeMesh():#(mesh_name, domain_dimensions, radius, scale):
#mesh_name = "3D_mesh"
    mesh_name = "many_holes_new"
    domain_dimensions = [1.0, 2.0]#, 1.0]
    radius = 0.2
    holes_x = 1
    holes_y = 2
    num_holes = holes_x * holes_y
    dim = len(domain_dimensions)
    scale = 1.0
    pythonToGeo(mesh_name, domain_dimensions, radius, holes_x, holes_y)
    #call(["gmsh", "-%i" % dim, "%s.geo" % mesh_name, "-clscale", "%f" % scale, "-format", "msh2", "-2"]);

#    return Mesh("%s.msh" % mesh_name, dim=dim)
    return (mesh_name, domain_dimensions, radius, num_holes)
 
