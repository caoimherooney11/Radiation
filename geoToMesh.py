from firedrake import Mesh
from subprocess import call
#from dolfin import mpi_comm_world, XDMFFile


def geoToMesh(code, name, scale, dim):
    with open("%s.geo" %name , "w") as text_file:
        text_file.write(code);

    call(["gmsh", "-%i" % dim, "%s.geo" % name, "-clscale", "%f" %scale, "-format", "msh2", "-2"]);
    return Mesh("%s.msh" % name, dim=dim)
