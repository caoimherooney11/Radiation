from firedrake import *
#from test_torusMesh import getTorusMesh
import numpy as np

def torusToSquare(mesh):
    Lx = 1
    Ly = 1
    
    coord_family = 'DG'
    coord_fs = VectorFunctionSpace(mesh, coord_family, 1, dim=2)
    old_coordinates = mesh.coordinates
    new_coordinates = Function(coord_fs)
    
    periodic_kernel = """
    double pi = 3.141592653589793;
    double eps = 1e-12;
    double bigeps = 1e-1;
    double phi, theta, Y, Z;
    Y = 0.0;
    Z = 0.0;
    
    for(int i=0; i<old_coords.dofs; i++) {
        Y += old_coords[i][1];
        Z += old_coords[i][2];
    }
    
    for(int i=0; i<new_coords.dofs; i++) {
        phi = atan2(old_coords[i][1], old_coords[i][0]);
        if (fabs(sin(phi)) > bigeps)
            theta = atan2(old_coords[i][2], old_coords[i][1]/sin(phi) - 1.0);
        else
            theta = atan2(old_coords[i][2], old_coords[i][0]/cos(phi) - 1.0);
    
        new_coords[i][0] = phi/(2.0*pi);
        if(new_coords[i][0] < -eps) {
            new_coords[i][0] += 1.0;
        }
        if(fabs(new_coords[i][0]) < eps && Y < 0.0) {
            new_coords[i][0] = 1.0;
        }
    
        new_coords[i][1] = theta/(2.0*pi);
        if(new_coords[i][1] < -eps) {
            new_coords[i][1] += 1.0;
        }
        if(fabs(new_coords[i][1]) < eps && Z < 0.0) {
            new_coords[i][1] = 1.0;
        }
    
        new_coords[i][0] *= Lx[0];
        new_coords[i][1] *= Ly[0];
    }
    """
    
    cLx = Constant(Lx)
    cLy = Constant(Ly)
    
    par_loop(periodic_kernel, dx,
             {"new_coords": (new_coordinates, WRITE),
              "old_coords": (old_coordinates, READ),
              "Lx": (cLx, READ),
              "Ly": (cLy, READ)})
    
    mesh_new = Mesh(new_coordinates)
    Vc = mesh_new.coordinates.function_space()
    x, y = SpatialCoordinate(mesh_new)
    f = Function(Vc).interpolate(as_vector([x-0.5, y-0.5]))
    mesh_new.coordinates.assign(f)

    File("test_mesh.pvd").write(mesh_new.coordinates)
#    mesh_new.Export("export.msh","Gmsh Format")
    return mesh_new

